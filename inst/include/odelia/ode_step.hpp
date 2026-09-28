// -*-c++-*-
#ifndef ODELIA_ODE_STEP_HPP_
#define ODELIA_ODE_STEP_HPP_

#include <array>
#include <vector>
#include <cstddef>
#include <XAD/XAD.hpp>
#include <odelia/adjoint.hpp>
#include <odelia/ode_interface.hpp>

namespace odelia {
namespace ode {

template <class System>
class Step {
public:
  using value_type = typename System::value_type;
  using state_type = std::vector<value_type>;
  
  // What one step's SIX rate evaluations solve for. One name, because the forward
  // walk, the sweep and the record all have to agree on the shape.
  //
  // ⚠️ SIX AND NOT FIVE, and the sixth is the one to understand. Five of them are
  // the stages; the sixth is the evaluation at the state the step ends at, which
  // first-same-as-last hands the next step as its own k1. A SWEEP re-derives that
  // one at the state it was handed and reads only 0..4 -- but a FORWARD replay
  // cannot re-derive it, because re-deriving is exactly what it is replaying to
  // avoid, and a step whose k1 was re-derived is wrong at first order in h.
  using solved_row = std::array<solved_values_t<System>, 6>;

  void resize(size_t size_);
  size_t order() const;
  // `solved` is the row this step is about to create: what its five stages solve
  // for goes in, and nothing is asked of the System about where it is.
  //
  // Or the row an earlier run already created, where a caller hands a CONST one:
  // the stages then LOAD what that run solved instead of solving again. Which of
  // the two happens is the constness of what was handed over and nothing else --
  // the rule `solved_scope` already follows one level down -- so there is no mode
  // here to keep, and none to hold a stale answer to.
  template <class Row>
  void step(System& system, Row& solved,
            double time, double step_size,
	    state_type &y,
	    state_type &yerr,
	    const state_type &dydt_in,
	    state_type &dydt_out);

  // The step transposed, for as many seeds as are handed in: one recording of the
  // whole step, swept once per seed. The active System is the walk's, held across
  // every step of one width. A caller wanting one row passes a batch of
  // one; there is no separate entry point for that, because a second signature
  // over the same recording is a second place for the seam between the state and
  // the parameter halves to be got wrong.
  void step_adjoint(active_system<System>& active,
                    const solved_row& solved,
                    double time, double step_size,
                    const state_type &y, const adjoint_rows& lambda_out,
                    adjoint_rows& lambda_in, adjoint_rows& parameter_adjoint);

  // Rate evaluations the sweeps since the last clear have recorded, counted where
  // they are recorded rather than added up as a total the loop could disagree
  // with. Six a step, whatever the seed count, because the step is recorded once
  // and swept per seed. A term entering once a step where it belongs once a stage
  // divides this by six, and no gradient check can see that, because a tangent and
  // a sweep apply the same multiplier.
  std::size_t recorded_rates = 0;

  static const bool can_use_dydt_in = true;
  static const bool first_same_as_last = true;

private:
  // The tableau, written once and used at whatever scalar the caller holds its
  // rates in: the forward step and the recording its transpose is taken from
  // both step through these, so the two cannot come apart.
  //
  // Y_i for stage i, into `out`: y at stage 0, and y plus the combination of the
  // earlier stage rates above that. Callers pass 1..5; see stage_row for why the
  // stage-0 arms stay.
  template <class S>
  void stage_state(int i, const std::vector<S>& y,
                   const std::vector<std::vector<S>>& k, double h,
                   std::vector<S>& out) const;
  // And the state the step ends at, y + h * (c1 k1 + c3 k3 + c4 k4 + c6 k6).
  // k2 and k5 reach it only through the later stages. `out` may be `y`.
  template <class S>
  void step_end(const std::vector<S>& y, const std::vector<std::vector<S>>& k,
                double h, std::vector<S>& out) const;
  double stage_time(int i, double time, double h) const;
  const double* stage_row(int i) const;
  // stage_state for a stage count known at compile time: the same sums in the
  // same order, with the earlier-stage count a constant so the inner loop
  // unrolls, and the data pointers held rather than re-read through the nested
  // vectors after every store. stage_state dispatches here for stages 1..5.
  template <int I, class S>
  void stage_state_fixed(const std::vector<S>& y,
                         const std::vector<std::vector<S>>& k, double h,
                         std::vector<S>& out) const;

  size_t size;
  std::vector<state_type> k{6};
  state_type ytmp;

  // Cash carp constants, from GSL.
  static const double ah[];
  static const double b21;
  static const double b3[];
  static const double b4[];
  static const double b5[];
  static const double b6[];
  static const double c1;
  static const double c3;
  static const double c4;
  static const double c6;

  // These are the differences of fifth and fourth order coefficients
  // for error estimation
  static const double ec[];
};

template <class System>
void Step<System>::resize(size_t size_) {
  size = size_;
  for (state_type& stage_rate : k) {
    stage_rate.resize(size);
  }
  ytmp.resize(size);
}

template <class System>
size_t Step<System>::order() const {
  // In GSL, comment says "FIXME: should this be 4?"
  return 5;
}

template <class System>
template <class Row>
void Step<System>::step(System& system,
                        Row& solved,
                        double time, double step_size,
                        state_type &y,
                        state_type &yerr,
                        const state_type &dydt_in,
                        state_type &dydt_out) {
  const double h = step_size;

  // First-same-as-last: k1 is the previous step's dydt_out, so the step costs five
  // rate evaluations and one more to hand the next step its own k1 -- which is why
  // that last one is addressed as the next step's stage 0.
  // A stage's rates, handed the slot it stores what it solves for into. A System
  // that solves for nothing is handed nothing and the branch compiles away.
  auto rates_at = [&](int i, const state_type& at, state_type& into) -> void {
    if constexpr (SolvesForValues<System>) {
      ode::derivs(system, at, into, stage_time(i, time, h), solved[i - 1]);
    } else {
      ode::derivs(system, at, into, stage_time(i, time, h));
    }
  };

  std::copy(dydt_in.begin(), dydt_in.end(), k[0].begin());
  // The stages written out, each at a compile-time stage count, so the kernel
  // inlines into this function: a run-time stage index kept it out of line, a
  // call per stage. The same kernels stage_state dispatches to, so the sweep
  // still reverses exactly this arithmetic.
  stage_state_fixed<1>(y, k, h, ytmp); rates_at(1, ytmp, k[1]);
  stage_state_fixed<2>(y, k, h, ytmp); rates_at(2, ytmp, k[2]);
  stage_state_fixed<3>(y, k, h, ytmp); rates_at(3, ytmp, k[3]);
  stage_state_fixed<4>(y, k, h, ytmp); rates_at(4, ytmp, k[4]);
  stage_state_fixed<5>(y, k, h, ytmp); rates_at(5, ytmp, k[5]);

  step_end(y, k, h, y);
  // The sixth evaluation, at the state the step ends at, which first-same-as-last
  // hands the next step as its own first rates -- so it carries a slot of its own:
  // see `solved_row` for why a sweep never reads it and a forward replay must.
  if constexpr (SolvesForValues<System>) {
    ode::derivs(system, y, dydt_out, time + h, solved[5]);
  } else {
    ode::derivs(system, y, dydt_out, time + h);
  }

  // Difference between 4th and 5th order, for error calculations
  const value_type* const k0 = k[0].data();
  const value_type* const k2 = k[2].data();
  const value_type* const k3 = k[3].data();
  const value_type* const k4 = k[4].data();
  const value_type* const k5 = k[5].data();
  for (size_t q = 0; q < size; ++q) {
    yerr[q] = h * (ec[1] * k0[q] + ec[3] * k2[q] + ec[4] * k3[q] +
                   ec[5] * k4[q] + ec[6] * k5[q]);
  }
}

// The tableau row stage i's state combines the earlier stage rates with.
//
// The stage-0 entry is kept although no caller passes 0, and so are the stage-0
// arms of stage_time and stage_state. Dropping them and re-basing this table at
// stage 2 reads as tidying away unreachable code, and it is not: this function is
// pure and inlined, so the compiler may evaluate it above stage_state's own
// `i == 1` early return, and `rows[i - 2]` is then an out-of-bounds read of a
// stack array at i == 1. It costs nothing to keep the index at i - 1 and one
// entry in the table, and the version that removed them returned a wrong
// gradient on the second of two calls.
template <class System>
const double* Step<System>::stage_row(int i) const {
  const double* const rows[] = {&b21, b3, b4, b5, b6};
  return rows[i - 1];
}

template <class System>
double Step<System>::stage_time(int i, double time, double h) const {
  return i == 0 ? time : time + ah[i - 1] * h;
}

template <class System>
template <class S>
void Step<System>::stage_state(int i, const std::vector<S>& y,
                               const std::vector<std::vector<S>>& k, double h,
                               std::vector<S>& out) const {
  if (i == 0) {
    std::copy(y.begin(), y.end(), out.begin());
    return;
  }
  switch (i) {
    case 1: return stage_state_fixed<1>(y, k, h, out);
    case 2: return stage_state_fixed<2>(y, k, h, out);
    case 3: return stage_state_fixed<3>(y, k, h, out);
    case 4: return stage_state_fixed<4>(y, k, h, out);
    default: return stage_state_fixed<5>(y, k, h, out);
  }
}

template <class System>
template <int I, class S>
void Step<System>::stage_state_fixed(const std::vector<S>& y,
                                     const std::vector<std::vector<S>>& k,
                                     double h, std::vector<S>& out) const {
  // Stage 1 keeps its single term grouped as b21 * h * k1: h * (b21 * k1)
  // rounds differently, and the reference numbers were blessed on this one.
  if constexpr (I == 1) {
    const S* const k0 = k[0].data();
    const S* const yp = y.data();
    S* const op = out.data();
    for (size_t q = 0; q < size; ++q) {
      op[q] = yp[q] + b21 * h * k0[q];
    }
    return;
  }
  const double* const b = stage_row(I);
  const S* kp[I];
  for (int m = 0; m < I; ++m) kp[m] = k[m].data();
  const S* const yp = y.data();
  S* const op = out.data();
  for (size_t q = 0; q < size; ++q) {
    // Summed in ascending stage, then one h. Cash-Karp's rows are dense, so
    // this is a sum over every earlier stage rather than a term for the
    // immediate predecessor.
    S combination = b[0] * kp[0][q];
    for (int m = 1; m < I; ++m) {
      combination += b[m] * kp[m][q];
    }
    op[q] = yp[q] + h * combination;
  }
}

template <class System>
template <class S>
void Step<System>::step_end(const std::vector<S>& y,
                            const std::vector<std::vector<S>>& k, double h,
                            std::vector<S>& out) const {
  const S* const k0 = k[0].data();
  const S* const k2 = k[2].data();
  const S* const k3 = k[3].data();
  const S* const k5 = k[5].data();
  const S* const yp = y.data();
  S* const op = out.data();
  for (size_t q = 0; q < size; ++q) {
    const S combination = c1 * k0[q] + c3 * k2[q] + c4 * k3[q] + c6 * k5[q];
    op[q] = yp[q] + h * combination;
  }
}

// lambda_in[m] = (d y_end / d y)^T lambda_out[m] for the one step step() takes
// from y, and the parameter rows alongside it.
//
// ONE recording spans the whole step: its six rate evaluations and the
// combination closing them. What the sweep transposes is therefore the
// arithmetic the stepper performs, where a recording per stage left the tableau
// to be transposed by hand beside the stepper and held consistent with it by
// discipline. The stage states are intermediates of the recording rather than a
// double rebuild ahead of it, so the step costs six model evaluations and not
// thirteen.
//
// The recording is derivs(), which is what the forward pass calls, so no System
// writes a transpose of its own; and the parameters ride in the same recording,
// so a stage the parameters reach carries their rows too.
template <class System>
void Step<System>::step_adjoint(active_system<System>& active,
                                const solved_row& solved,
                                double time, double step_size,
                                const state_type &y, const adjoint_rows& lambda_out,
                                adjoint_rows& lambda_in,
                                adjoint_rows& parameter_adjoint) {
  using scalar = active_scalar<double>;
  const double h = step_size;
  if (lambda_out.empty()) {
    util::stop("step_adjoint: needs at least one seed");
  }
  // The recording hands the whole state buffer to the slice below, and `size` is
  // what resize() set -- so a state of another width is checked here rather than
  // in the two callers above that happen to check it.
  util::check_length(y.size(), size);

  auto whole_step = [&](auto& sys,
                        typename std::vector<scalar>::const_iterator x,
                        std::vector<scalar>& y_end) -> void {
    const std::vector<scalar> y0(x, x + static_cast<std::ptrdiff_t>(size));
    std::vector<std::vector<scalar>> rate(6, std::vector<scalar>(size));
    std::vector<scalar> stage(size);
    // k1 is re-derived at this step's own start state, and unaddressed on purpose:
    // the run took its first rates either at the end of the step before this one or,
    // where it widened in between, at a state no record holds. A descent that starts
    // at an arbitrary step cannot tell those apart, so it asks for neither.
    ode::derivs(sys, y0, rate[0], time);
    ++recorded_rates;
    for (int i = 1; i < 6; ++i) {
      stage_state(i, y0, rate, h, stage);
      if constexpr (SolvesForValues<System>) {
        ode::derivs(sys, stage, rate[i], stage_time(i, time, h),
                    std::as_const(solved[i - 1]));
      } else {
        ode::derivs(sys, stage, rate[i], stage_time(i, time, h));
      }
      ++recorded_rates;
    }
    step_end(y0, rate, h, y_end);
  };

  ode::state_and_parameter_adjoints(active, y, lambda_out, whole_step, lambda_in,
                                    parameter_adjoint);
}

// RKCK coefficients, from GSL
template <class System>
const double Step<System>::ah[] = {
  1.0 / 5.0, 0.3, 3.0 / 5.0, 1.0, 7.0 / 8.0 };

template <class System>
const double Step<System>::b21 = 1.0 / 5.0;
template <class System>
const double Step<System>::b3[] = { 3.0 / 40.0, 9.0 / 40.0 };
template <class System>
const double Step<System>::b4[] = { 0.3, -0.9, 1.2 };
template <class System>
const double Step<System>::b5[] = {
  -11.0 / 54.0, 2.5, -70.0 / 27.0, 35.0 / 27.0 };

template <class System>
const double Step<System>::b6[] = {
  1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0,
  44275.0 / 110592.0, 253.0 / 4096.0 };

template <class System>
const double Step<System>::c1 = 37.0 / 378.0;
template <class System>
const double Step<System>::c3 = 250.0 / 621.0;
template <class System>
const double Step<System>::c4 = 125.0 / 594.0;
template <class System>
const double Step<System>::c6 = 512.0 / 1771.0;

template <class System>
const double Step<System>::ec[] = {
  0.0, 37.0 / 378.0 - 2825.0 / 27648.0, 0.0,
  250.0 / 621.0 - 18575.0 / 48384.0,
  125.0 / 594.0 - 13525.0 / 55296.0,
  -277.0 / 14336.0, 512.0 / 1771.0 - 0.25 };

}
}

#endif
