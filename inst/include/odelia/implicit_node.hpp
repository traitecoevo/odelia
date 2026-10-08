// -*-c++-*-
#ifndef ODELIA_IMPLICIT_NODE_HPP_
#define ODELIA_IMPLICIT_NODE_HPP_

#include <cmath>
#include <cstddef>
#include <span>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <XAD/XAD.hpp>
#include <odelia/ode_interface.hpp>
#include <odelia/ode_util.hpp>

namespace odelia {

// One input a value responds to, and how. Kept as a pair so the two cannot be
// assembled from separate lists and paired by position.
//
// THE INPUT IS BORROWED, and must outlive the record it is handed to. Held by
// reference because copying an active scalar is not free: the copy registers a
// tape slot and records a statement, for a value that is only ever read: n
// inputs would cost n of those per row.
template <class S>
struct input_and_derivative {
  const S& input;
  double derivative;
};

// What a record could not record, beside the value it wrote. A row that cannot be
// recorded is not an error here: the value is still the value, and whether a
// consumer can go on without the row is the consumer's to decide -- one output's
// rows can go missing while another's survive, and a stop takes both.
struct record_report {
  bool whole = true;
  // Which input's row is missing, and what was wrong with it. Meaningless where
  // `whole`.
  std::size_t at = 0;
  std::string why;
};

// The scalar that carries an adjoint, so a row set can be a statement's operand
// run. A direction has no tape to hold one and takes the arithmetic below.
template <class S>
concept CarriesAdjoint = xad::ExprTraits<S>::isReverse;

// Cost, in tape statements walked per solve of a submodel of T statements with
// m outputs, read by a consumer sweeping k seeds:
//
//   recorded inline on the consumer's tape   T + kT
//   supplied here, nothing recorded          m + km
//
// The second carries no T, so what a supplied row costs does not depend on how
// long the submodel is. Recording the submodel on a tape of its own and sweeping
// it m times to extract a dense block pays only when the submodel has fewer
// outputs than the consumer has seeds.
//
// `into` receives `value` carrying the derivatives supplied against it: the
// number is `value` itself, and its derivative with respect to each input is the
// one supplied.
//
// This is how a quantity computed away from the tape gets onto it -- a
// root-find, a submodel's own solve, anything whose derivative is known by some
// means other than recording the steps that produced it. At a plain double the
// rows vanish and this is the value.
//
// ONE STATEMENT, whatever the row count. A tape statement is a left-hand side
// over a run of operations, so n rows are n operations under one lhs -- not the
// n recorded assignments that writing the sum out as `out += d * (x -
// to_passive(x))` costs.
//
// NOTHING PARTIAL. Every row is tested before any is recorded, because a value
// carrying some of its rows is a channel that has gone missing with every number
// still finite -- which is worse than no rows at all, since a consumer told the
// rows are absent can carry the value as a constant and say so.
//
// The report is the return, so a caller cannot take the number without being
// handed the reading of it.
template <class S>
[[nodiscard]] record_report record_with_derivatives(
    double value, std::span<const input_and_derivative<S>> against, S& into) {
  for (std::size_t i = 0; i < against.size(); ++i) {
    // Either half being non-finite poisons the VALUE and not only what is
    // recorded against it: NaN times zero is not a number, and an infinite input
    // minus its own passive copy is NaN rather than zero. Tested here, where the
    // halves are still separable; downstream they are one expression.
    const double at = util::to_passive(against[i].input);
    if (!util::is_finite(against[i].derivative) || !util::is_finite(at)) {
      into = value;
      return {false, i,
              "input " + util::to_string(static_cast<int>(i)) + " of " +
                  util::to_string(static_cast<int>(against.size())) +
                  " has value " + util::format_double(at) + " and derivative " +
                  util::format_double(against[i].derivative) +
                  ", one of which is not finite, so the value they belong to "
                  "cannot be recorded"};
    }
  }
  // A zero derivative contributes exactly zero to the value and exactly nothing to
  // the transpose, and carrying it costs a tape edge the sweep then walks.
  if constexpr (CarriesAdjoint<S>) {
    using tape_type = typename S::tape_type;
    // ⚠️ THE ORDER IS THE VALUE, THEN THE ROWS, THEN THE CLOSE, and each step is
    // load-bearing. A statement's operations are everything pushed since the last
    // left-hand side, so anything that pushes between the first row and the close
    // claims them -- which is why every row is tested above rather than here.
    // `registerOutput` is a no-op on a value that already holds a slot, leaving
    // the rows for whatever statement closes next, so the destination is built
    // here and moved out.
    //
    // Pushed one row at a time rather than gathered first: operations accumulate
    // until the left-hand side closes over them, so the run is the same run and
    // the two scratch arrays a gather needs are two allocations per call at a
    // rate of millions.
    S out = value;
    if (tape_type* tape = tape_type::getActive()) {
      for (const input_and_derivative<S>& term : against) {
        if (term.derivative == 0.0) {
          continue;
        }
        // A passive input holds no slot, and the sweep indexes the slot it is
        // pushed without a bounds check, so pushing one would corrupt memory
        // rather than raise. It has no row to carry either way.
        const typename tape_type::slot_type slot = term.input.getSlot();
        if (slot == tape_type::INVALID_SLOT) {
          continue;
        }
        tape->pushAll(&term.derivative, &slot, 1u);
      }
      tape->registerOutput(out);
    }
    into = std::move(out);
  } else {
    S out = value;
    for (const input_and_derivative<S>& term : against) {
      if (term.derivative == 0.0) {
        continue;
      }
      // Zero in value and carrying the derivative, so the number is untouched.
      out += term.derivative * (term.input - util::to_passive(term.input));
    }
    into = out;
  }
  return {};
}
// The value y* defined implicitly by a scalar equation F(y) = 0, made
// differentiable. y* is solved in double, off the tape, by whatever root-find the
// caller already has; this returns it on the tape carrying the derivative the
// implicit function theorem gives:
//     dy*/dp = -(dF/dp) / (dF/dy).
//
// dF/dp is taped: F is evaluated once, at S, and its derivative in every active
// input the residual reads IS dF/dp. dF/dy is the caller's, because a residual
// mostly knows its own slope: it is
// one term of the expression the caller just wrote. One that does not can take a
// tangent through its own body and hand the answer here, which is a local three
// lines rather than a rebindable parameter object every caller must carry.
//
// dF/dy is read at the point, so its own dependence on p belongs to the second
// derivative rather than to this one. It is what the theorem divides by, so a fold
// -- where it approaches zero and the quotient is garbage rather than large -- is
// the one thing this stops on.
template <class S, class Residual>
S implicit_value(double y_star, double dFdy, Residual&& F) {
  if constexpr (std::is_same_v<S, double>) {
    return y_star;
  } else {
    // A residual written `[](auto y) { return ...; }` deduces an XAD
    // expression-template return type holding references to the temporaries of
    // its return statement; those die when it returns, so evaluating it here
    // would read a destroyed tape slot -- a segfault on the reverse sweep, not a
    // wrong number. Requiring the scalar back exactly makes that a compile error.
    static_assert(
        std::is_same_v<std::invoke_result_t<Residual&, const S&>, S>,
        "implicit_value: the residual must return its own scalar exactly "
        "([](const S& y) -> S { ... }). A deduced return type is an expression "
        "template referencing temporaries that are dead by the time this "
        "evaluates it.");
    if (!util::is_finite(dFdy) || dFdy == 0.0) {
      util::stop("implicit_value: dF/dy is " + util::format_double(dFdy) +
                 " at the operating point, so the implicit function theorem "
                 "does not apply there (a fold?)");
    }
    // corr's value is ~0, since y* is the root; its derivative is (dF/dp)/(dF/dy),
    // so y* against it with a coefficient of -1 is the theorem's own quotient.
    std::vector<input_and_derivative<S>> corrections;
    const S corr = F(S(y_star)) / dFdy;
    corrections.push_back({corr, -1.0});
    S out;
    const record_report report =
        record_with_derivatives<S>(y_star, corrections, out);
    if (!report.whole) {
      util::stop("implicit_value: " + report.why);
    }
    return out;
  }
}


// The same value, with the residual's own statements taken OFF the caller's tape.
//
// The three-argument form above leaves them there: `corr` has to stay reachable
// for the outer sweep to walk it, so a residual costing T statements costs the
// consumer T statements per solve, and the consumer sweeps them once per output
// it asks for. This form sweeps the residual ONCE, here, keeps the numbers that
// fall out, and rewinds. The consumer sees one statement.
//
// `inputs` are the values whose rows are wanted, in whatever shape they are held
// -- `visit_active` opens a scalar, a container, a pointer, or anything declaring
// `for_each_active`. They are the residual's own inputs, so a caller hands over
// what it already holds rather than assembling a list.
//
// ⚠️ A SHAPE `visit_active` DOES NOT OPEN IS SKIPPED IN SILENCE (see visit_active),
// AND THE COLUMN IT COSTS DOES NOT COME BACK ZERO. The sweep below deposits the
// missing row on that input's slot, and nothing here can reach it to clear it:
// XAD's `clearDerivativesAfter` reaches only slots created after the mark, and an
// input predates it by construction. So the consumer's own sweep adds to what was
// left. Nothing here can detect it either, so no count is handed back: see
// `ode::count_active_slots` for what a test that knows its own list can assert.
//
// The caller's own adjoints on the inputs are left as they were.
//
// At least one input, so the arity says what the form needs: rewinding with none
// would discard the residual and supply nothing in its place. Three arguments is
// the other overload, which leaves the residual on the tape.
template <class S, class Residual, class First, class... Rest>
S implicit_value(double y_star, double dFdy, Residual&& F, const First& first,
                 const Rest&... rest) {
  // A direction has no tape to rewind, and at a double there is nothing to
  // record at all: both take the arithmetic form, which carries the same rows
  // where it has any.
  if constexpr (!CarriesAdjoint<S>) {
    return implicit_value<S>(y_star, dFdy, std::forward<Residual>(F));
  } else {
    static_assert(
        std::is_same_v<std::invoke_result_t<Residual&, const S&>, S>,
        "implicit_value: the residual must return its own scalar exactly "
        "([](const S& y) -> S { ... }). A deduced return type is an expression "
        "template referencing temporaries that are dead by the time this "
        "evaluates it.");
    if (!util::is_finite(dFdy) || dFdy == 0.0) {
      util::stop("implicit_value: dF/dy is " + util::format_double(dFdy) +
                 " at the operating point, so the implicit function theorem "
                 "does not apply there (a fold?)");
    }
    using tape_type = typename S::tape_type;
    tape_type* tape = tape_type::getActive();
    if (tape == nullptr) {
      return S(y_star);
    }
    const typename tape_type::position_type mark = tape->getPosition();

    // ⚠️ THE INPUTS' ADJOINTS ARE HELD AND PUT BACK, NOT ZEROED. The sweep below
    // runs on the CALLER'S tape and ACCUMULATES, so a slot the caller has already
    // swept something into arrives non-zero -- reading it afterwards would report
    // the caller's own answer as this residual's row, and zeroing it afterwards
    // would destroy the caller's.
    //
    // Gathered before the sweep for the same reason: the slot list has to be the
    // one that was zeroed, not one read back out of a tape the sweep has touched.
    std::vector<typename tape_type::slot_type> slots;
    auto gather = [&](const S& x) {
      const typename tape_type::slot_type slot = x.getSlot();
      if (slot != tape_type::INVALID_SLOT) {
        slots.push_back(slot);
      }
    };
    odelia::ode::visit_active(gather, first, rest...);
    std::vector<double> held(slots.size());
    for (std::size_t i = 0; i < slots.size(); ++i) {
      held[i] = tape->derivative(slots[i]);
      tape->derivative(slots[i]) = 0.0;
    }

    {
      S corr = F(S(y_star)) / dFdy;
      tape->registerOutput(corr);
      xad::derivative(corr) = 1.0;
      tape->computeAdjointsTo(mark);
    }
    // Rewound BEFORE the rows are written, so the statement closed below carries
    // them rather than the residual's own.
    tape->resetTo(mark);
    S out = y_star;
    // Read and cleared through the TAPE, by slot, so an input can arrive const --
    // which is how every caller already holds the things a residual reads.
    for (std::size_t i = 0; i < held.size(); ++i) {
      // What the slot holds now is the residual's deposit alone; the caller's
      // own adjoint goes back.
      const double adj = tape->derivative(slots[i]);
      tape->derivative(slots[i]) = held[i];
      const double row = -adj;
      if (row == 0.0) {
        continue;
      }
      // -(dF/dp)/(dF/dy) is the theorem; dF/dy divided the residual above.
      tape->pushAll(&row, &slots[i], 1u);
    }
    tape->registerOutput(out);
    return out;
  }
}


}

#endif