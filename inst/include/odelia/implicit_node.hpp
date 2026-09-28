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
// to_passive(x))` costs. Carrying a leaf-shaped boundary that way put the whole
// of a submodel's arithmetic on a consumer's tape.
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
        // ⚠️ A PASSIVE INPUT HOLDS NO SLOT, and the sweep indexes the slot it is
        // pushed without a bounds check, so pushing one corrupts memory rather
        // than raising. It has no row to carry either way: nothing outside reads
        // it.
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
//
// ⚠️ THIS ONCE HAD A REPORTING SIBLING, and the sibling went because the distinction
// was one nothing acted on. The argument for it was that a BOUND must stop -- its
// value IS what the equation defines -- while an INTERIOR optimum could carry on with
// the row missing, since the envelope theorem spares the objective. Both were true of
// the theorem and neither was true of the consumer: the report's one reader turned it
// straight into the same metric-level refusal the catch around this makes. Two
// mechanisms, one outcome, and a record_report threaded through LeafOutputs and two
// signatures to carry the difference.
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
// ⚠️ A SHAPE `visit_active` DOES NOT OPEN IS SKIPPED IN SILENCE, AND THE COLUMN
// IT COSTS DOES NOT COME BACK ZERO. The sweep below deposits the missing row on
// that input's slot, and nothing here can reach it to clear it: XAD's
// `clearDerivativesAfter` reaches only slots created after the mark, and an
// input predates it by construction. So the consumer's own sweep adds to what
// was left. Measured on a region whose true dw/dy is 27, with y undeclared: the
// answer is 11, not 7 and not 0.
//
// ⚠️ AND NOTHING HERE CAN DETECT IT, which is why this does not hand back a count
// for a caller to check. The only count available is of the values the list DOES
// name, so it says nothing about one left out -- and it cannot say a row was
// written either, since a named input whose adjoint comes back zero is still one
// of them. A test that knows its own list and wants the number calls
// `ode::count_active_slots` beside the call, where it reads as the assertion it
// is; five call sites carrying an out-parameter none of them reads does not.
//
// ⚠️ THE INPUTS' ADJOINTS ARE HELD AND PUT BACK, not zeroed, and the row is the
// difference. They accumulate on the caller's tape, so a second solve against the
// same inputs would otherwise add to the first -- and a caller that had already
// swept something into one of those slots would see its own answer reported as
// this node's row, and then destroyed.
// ⚠️ AT LEAST ONE INPUT, SPELLED IN THE SIGNATURE. Rewinding with none discards
// the residual and supplies nothing in its place, which is a value with no rows
// rather than a value -- so the arity says what the form needs. Three arguments
// is the other overload, which leaves the residual on the tape.
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
    // swept something into arrives non-zero -- and reading it afterwards reports
    // the caller's own answer as this residual's row. Measured on a region whose
    // true row is 12, against a caller that had left 9 on the same slot: the
    // harvest read 21, and zeroing afterwards destroyed the 9 as well.
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
      const double adj = tape->derivative(slots[i]);
      // The caller's own adjoint goes back, because this node had no business
      // taking it: what the residual deposited is adj - held[i].
      tape->derivative(slots[i]) = held[i];
      const double row = -(adj - held[i]);
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


// A whole region of a recording, replaced on the caller's tape by the rows of
// what it produced.
//
// `implicit_value` above does this for a region with ONE output, where the row
// set is the implicit function theorem's. This is the same trade with the outputs
// counted: the region is swept once per output rather than left for the consumer
// to sweep once per output it asks for, and what the consumer keeps is m
// statements instead of the region's whole arithmetic. It pays whenever the
// consumer sweeps more often than the region has outputs.
//
// `body()` runs the region and RETURNS the actives it produced -- returns rather
// than fills, because a region sizes its own outputs and a span taken before it
// runs points at whatever the sizing moved. `inputs` are the values whose rows
// are wanted, opened by `visit_active`.
//
// ⚠️ THE INPUTS MUST BE AN ANTICHAIN -- none computed from another. A row is
// attached for each of them, and the consumer's own sweep then carries each row
// onward through whatever produced it, so an input reachable from another input
// is counted once directly and again through its parent. Every number stays
// finite and the rows are quietly too large. The safe cut is the values as they
// crossed into the region, which is why `body` starts after they are all built.
//
// ⚠️ A SHAPE `visit_active` DOES NOT OPEN IS SKIPPED IN SILENCE -- a row that
// never arrives. No count comes back to say so, for implicit_value's reason: the
// only count available is of the values the list DOES name. A test that wants it
// calls `ode::count_active_slots` beside this.
//
// `scratch` is the caller's and is reused: this runs per cohort per stage per
// step, so an allocation here is an allocation there.
template <class S, class Body, class... Inputs>
void preaccumulate(Body&& body, std::vector<double>& scratch,
                   const Inputs&... inputs) {
  // A direction has no tape to rewind and carries its rows in the arithmetic, so
  // there the region simply runs.
  if constexpr (!CarriesAdjoint<S>) {
    body();
    return;
  } else {
  using tape_type = typename S::tape_type;
  tape_type* tape = tape_type::getActive();
  if (tape == nullptr) {
    // The scalar carries an adjoint, so the caller means to record: with no tape
    // running there is nothing to take the region off, and carrying on would
    // hand back a value whose rows are all missing.
    util::stop("preaccumulate: the scalar carries an adjoint but no tape is "
               "active, so there is no recording to take this region out of");
  }

  const typename tape_type::position_type mark = tape->getPosition();
  const std::span<S* const> outputs = body();

  // The slots to harvest, gathered once: the walk is over the caller's shapes
  // and costs more than reading a vector back.
  std::vector<typename tape_type::slot_type> slots;
  slots.clear();
  auto gather = [&](const S& x) {
    const typename tape_type::slot_type slot = x.getSlot();
    if (slot != tape_type::INVALID_SLOT) {
      slots.push_back(slot);
    }
  };
  odelia::ode::visit_active(gather, inputs...);

  // Held and put back, for the reason implicit_value gives: the sweeps below run
  // on the CALLER'S tape and accumulate, so a slot it has already swept into
  // arrives non-zero and would be read as this region's row.
  std::vector<double> held(slots.size());
  for (std::size_t i = 0; i < slots.size(); ++i) {
    held[i] = tape->derivative(slots[i]);
    tape->derivative(slots[i]) = 0.0;
  }

  const std::size_t m = outputs.size();
  const std::size_t n = slots.size();
  scratch.assign(m * n, 0.0);
  std::vector<double> values(m, 0.0);
  for (std::size_t j = 0; j < m; ++j) {
    values[j] = util::to_passive(*outputs[j]);
    tape->registerOutput(*outputs[j]);
    xad::derivative(*outputs[j]) = 1.0;
    tape->computeAdjointsTo(mark);
    // ⚠️ THE SEED IS CLEARED TOO. `resetTo` rewinds the statements and NOT the
    // slot counter, so a seeded output leaves a slot behind carrying 1.0; a later
    // sweep of the same recording reads it and answers differently depending on
    // where it started, which is what a decomposition check sees and a single
    // sweep does not.
    xad::derivative(*outputs[j]) = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      scratch[j * n + i] = tape->derivative(slots[i]);
      // Zeroed as it is read: the adjoints accumulate, so the next output would
      // otherwise carry this one's rows as well as its own.
      tape->derivative(slots[i]) = 0.0;
    }
  }

  tape->resetTo(mark);
  for (std::size_t i = 0; i < n; ++i) {
    tape->derivative(slots[i]) = held[i];
  }
  for (std::size_t j = 0; j < m; ++j) {
    S out = values[j];
    for (std::size_t i = 0; i < n; ++i) {
      const double row = scratch[j * n + i];
      if (row == 0.0) {
        continue;
      }
      tape->pushAll(&row, &slots[i], 1u);
    }
    tape->registerOutput(out);
    *outputs[j] = std::move(out);
  }
  }
}

}

#endif