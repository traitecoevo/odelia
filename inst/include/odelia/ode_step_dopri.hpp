// -*-c++-*-
#ifndef ODELIA_ODE_STEP_DOPRI_HPP_
#define ODELIA_ODE_STEP_DOPRI_HPP_

// Dormand-Prince 5(4): a 7-stage explicit Runge-Kutta pair of order 5 with an
// embedded 4th-order error estimate, first-same-as-last (the seventh stage is
// the derivative at the new point, which is the first stage of the next step),
// so six evaluations per step. The method behind MATLAB's ode45, deSolve's
// ode45 (rk45dp7), scipy's RK45 and Hairer's dopri5.f. Coefficients from
// Dormand & Prince (1980) as tabulated in Hairer, Norsett & Wanner, Solving
// ODEs I, table 5.2.
//
// What it has that the Cash-Karp pair (ode_step_rkck.hpp) has not is a free
// continuous extension of order 4 (Hairer's dopri5.f, `contd5`): the state
// anywhere inside an accepted step, from the step's own stages, at no further
// evaluation and with an error of the same order as the step's. That is what
// makes dense output honest -- output at times the stepper did not land on,
// still within the tolerance it was asked for -- where a cubic Hermite on the
// endpoints is one order short. See SolverInternal::interpolate and
// Solver::advance_collect.
//
// Otherwise it is the Cash-Karp stepper's interface member for member: a step
// writes or reads a solved_row (stages 2..6, then the seventh evaluation at the
// state the step ends at), and step_adjoint transposes one step, so a run under
// either is recorded, replayed and swept the same way. `dense()` is the extra.

#include <array>
#include <cstddef>
#include <utility>
#include <vector>
#include <XAD/XAD.hpp>
#include <odelia/adjoint.hpp>
#include <odelia/ode_interface.hpp>

namespace odelia {
namespace ode {

template <class System>
class DopriStep {
public:
  using value_type = typename System::value_type;
  using state_type = std::vector<value_type>;

  // Sized, not cleared: a sweep at the run's own width resizes on the way in,
  // and the last accepted step's stages are what dense() reads afterwards.
  void resize(size_t size_) {
    size = size_;
    for (state_type& stage_rate : k) {
      stage_rate.resize(size);
    }
    ytmp.resize(size);
    y0.resize(size);
  }

  // Order of the propagated solution (the controller uses 1/order).
  size_t order() const { return 5; }

  // `solved` is the row this step creates, or, where a caller hands a CONST one,
  // the row an earlier run created, which the stages then load: the same rule as
  // Step::step.
  template <class Row>
  void step(System& system, Row& solved, double time, double step_size,
            state_type& y, state_type& yerr,
            const state_type& dydt_in, state_type& dydt_out) {
    const double h = step_size;
    y0 = y;
    h_last = h;
    k[0] = dydt_in;

    stages(y, k, h, ytmp,
           [&](int i, const state_type& at, state_type& into) -> void {
             if constexpr (SolvesForValues<System>) {
               ode::derivs(system, at, into, stage_time(i, time, h),
                           solved[i - 1]);
             } else {
               ode::derivs(system, at, into, stage_time(i, time, h));
             }
           });

    // The fifth-order solution; its weights are the seventh stage's row.
    step_end(y0, k, h, y);
    // The seventh evaluation, in the row's sixth slot: see solved_row.
    if constexpr (SolvesForValues<System>) {
      ode::derivs(system, y, k[6], time + h, solved[5]);
    } else {
      ode::derivs(system, y, k[6], time + h);
    }
    dydt_out = k[6];

    // Difference between the fifth- and fourth-order solutions.
    for (size_t i = 0; i < size; ++i) {
      yerr[i] = h * (e1 * k[0][i] + e3 * k[2][i] + e4 * k[3][i] +
                     e5 * k[4][i] + e6 * k[5][i] + e7 * k[6][i]);
    }
  }

  // The step transposed, for as many seeds as are handed in; Step::step_adjoint
  // says what the arguments are. The seventh evaluation is not recorded: the
  // state the step ends at does not depend on it, only the error estimate and
  // the dense output do.
  void step_adjoint(active_system<System>& active,
                    const solved_row<System>& solved,
                    const carried_rates<System>* k1,
                    double time, double step_size,
                    const state_type& y, const adjoint_rows& lambda_out,
                    adjoint_rows& lambda_in, adjoint_rows& parameter_adjoint);
  void step_adjoint(active_system<System>& active,
                    const solved_row<System>& solved,
                    double time, double step_size, const state_type& y,
                    const adjoint_rows& lambda_out, adjoint_rows& lambda_in,
                    adjoint_rows& parameter_adjoint) {
    step_adjoint(active, solved, nullptr, time, step_size, y, lambda_out,
                 lambda_in, parameter_adjoint);
  }

  // Rate evaluations the sweeps since the last clear have recorded; see
  // Step::recorded_rates. Six a step here too.
  std::size_t recorded_rates = 0;

  // The state at fraction theta in [0, 1] of the last step taken, from that
  // step's stages (dopri5.f, contd5). Valid once a step has been taken at
  // this size; the caller knows whether that step was accepted.
  void dense(double theta, const state_type& y1, state_type& out) const {
    out.resize(size);
    const double h = h_last;
    const double t1 = 1.0 - theta;
    for (size_t i = 0; i < size; ++i) {
      const value_type r1 = y0[i];
      const value_type r2 = y1[i] - y0[i];
      const value_type r3 = h * k[0][i] - r2;
      const value_type r4 = r2 - h * k[6][i] - r3;
      const value_type r5 = h * (d1 * k[0][i] + d3 * k[2][i] + d4 * k[3][i] +
                                 d5 * k[4][i] + d6 * k[5][i] + d7 * k[6][i]);
      out[i] = r1 + theta * (r2 + t1 * (r3 + theta * (r4 + t1 * r5)));
    }
  }

  static const bool can_use_dydt_in = true;
  static const bool first_same_as_last = true;

private:
  // The tableau, written once and used at whatever scalar the caller holds its
  // rates in: the forward step and the recording its transpose is taken from
  // both step through these, so the two agree to the last place where the
  // compiler does not contract (see Step's kernels).
  //
  // Stages 2..6 from y and k[0], each stage's rates into k[i] by `rates(i, at,
  // into)`. One walk for both passes, so the unrolled stages exist once.
  template <class S, class Rates>
  void stages(const std::vector<S>& y, std::vector<std::vector<S>>& k_,
              double h, std::vector<S>& tmp, Rates&& rates) const {
    stage_state<1>(y, k_, h, tmp); rates(1, tmp, k_[1]);
    stage_state<2>(y, k_, h, tmp); rates(2, tmp, k_[2]);
    stage_state<3>(y, k_, h, tmp); rates(3, tmp, k_[3]);
    stage_state<4>(y, k_, h, tmp); rates(4, tmp, k_[4]);
    stage_state<5>(y, k_, h, tmp); rates(5, tmp, k_[5]);
  }

  // Y for stage I + 1, into `out`. Each stage keeps its own written-out sum,
  // the one the stepper's results were taken with: a running sum is fused
  // differently under contraction and would move them.
  template <int I, class S>
  void stage_state(const std::vector<S>& y,
                   const std::vector<std::vector<S>>& k_, double h,
                   std::vector<S>& out) const {
    const S* const yp = y.data();
    const S* const k0 = k_[0].data();
    const S* const k1 = k_[1].data();
    const S* const k2 = k_[2].data();
    const S* const k3 = k_[3].data();
    const S* const k4 = k_[4].data();
    S* const op = out.data();
    for (size_t q = 0; q < size; ++q) {
      if constexpr (I == 1) {
        op[q] = yp[q] + h * a21 * k0[q];
      } else if constexpr (I == 2) {
        op[q] = yp[q] + h * (a31 * k0[q] + a32 * k1[q]);
      } else if constexpr (I == 3) {
        op[q] = yp[q] + h * (a41 * k0[q] + a42 * k1[q] + a43 * k2[q]);
      } else if constexpr (I == 4) {
        op[q] = yp[q] + h * (a51 * k0[q] + a52 * k1[q] + a53 * k2[q] +
                             a54 * k3[q]);
      } else {
        op[q] = yp[q] + h * (a61 * k0[q] + a62 * k1[q] + a63 * k2[q] +
                             a64 * k3[q] + a65 * k4[q]);
      }
    }
  }

  // The state the step ends at. k[1] reaches it only through the later stages.
  // `out` may be `y`.
  template <class S>
  void step_end(const std::vector<S>& y,
                const std::vector<std::vector<S>>& k_, double h,
                std::vector<S>& out) const {
    const S* const yp = y.data();
    const S* const k0 = k_[0].data();
    const S* const k2 = k_[2].data();
    const S* const k3 = k_[3].data();
    const S* const k4 = k_[4].data();
    const S* const k5 = k_[5].data();
    S* const op = out.data();
    for (size_t q = 0; q < size; ++q) {
      op[q] = yp[q] + h * (b1 * k0[q] + b3 * k2[q] + b4 * k3[q] +
                           b5 * k4[q] + b6 * k5[q]);
    }
  }

  double stage_time(int i, double time, double h) const {
    switch (i) {
      case 1: return time + c2 * h;
      case 2: return time + c3 * h;
      case 3: return time + c4 * h;
      case 4: return time + c5 * h;
      default: return time + h;
    }
  }

  size_t size = 0;
  std::vector<state_type> k{7};
  state_type ytmp;
  state_type y0;
  double h_last = 0.0;

  static constexpr double c2 = 1.0 / 5.0;
  static constexpr double c3 = 3.0 / 10.0;
  static constexpr double c4 = 4.0 / 5.0;
  static constexpr double c5 = 8.0 / 9.0;
  static constexpr double a21 = 1.0 / 5.0;
  static constexpr double a31 = 3.0 / 40.0;
  static constexpr double a32 = 9.0 / 40.0;
  static constexpr double a41 = 44.0 / 45.0;
  static constexpr double a42 = -56.0 / 15.0;
  static constexpr double a43 = 32.0 / 9.0;
  static constexpr double a51 = 19372.0 / 6561.0;
  static constexpr double a52 = -25360.0 / 2187.0;
  static constexpr double a53 = 64448.0 / 6561.0;
  static constexpr double a54 = -212.0 / 729.0;
  static constexpr double a61 = 9017.0 / 3168.0;
  static constexpr double a62 = -355.0 / 33.0;
  static constexpr double a63 = 46732.0 / 5247.0;
  static constexpr double a64 = 49.0 / 176.0;
  static constexpr double a65 = -5103.0 / 18656.0;
  static constexpr double b1 = 35.0 / 384.0;
  static constexpr double b3 = 500.0 / 1113.0;
  static constexpr double b4 = 125.0 / 192.0;
  static constexpr double b5 = -2187.0 / 6784.0;
  static constexpr double b6 = 11.0 / 84.0;
  // b - bhat: the error estimate's weights.
  static constexpr double e1 = 71.0 / 57600.0;
  static constexpr double e3 = -71.0 / 16695.0;
  static constexpr double e4 = 71.0 / 1920.0;
  static constexpr double e5 = -17253.0 / 339200.0;
  static constexpr double e6 = 22.0 / 525.0;
  static constexpr double e7 = -1.0 / 40.0;
  // Dense output (dopri5.f).
  static constexpr double d1 = -12715105075.0 / 11282082432.0;
  static constexpr double d3 = 87487479700.0 / 32700410799.0;
  static constexpr double d4 = -10690763975.0 / 1880347072.0;
  static constexpr double d5 = 701980252875.0 / 199316789632.0;
  static constexpr double d6 = -1453857185.0 / 822651844.0;
  static constexpr double d7 = 69997945.0 / 29380423.0;
};

// lambda_in[m] = (d y_end / d y)^T lambda_out[m] for the one step step() takes
// from y, and the parameter rows alongside it: one recording of the step's six
// rate evaluations and the combination closing them, as Step::step_adjoint.
template <class System>
void DopriStep<System>::step_adjoint(active_system<System>& active,
                                     const solved_row<System>& solved,
                                     const carried_rates<System>* k1,
                                     double time, double step_size,
                                     const state_type& y,
                                     const adjoint_rows& lambda_out,
                                     adjoint_rows& lambda_in,
                                     adjoint_rows& parameter_adjoint) {
  using scalar = active_scalar<double>;
  const double h = step_size;
  if (lambda_out.empty()) {
    util::stop("step_adjoint: needs at least one seed");
  }
  util::check_length(y.size(), size);

  auto whole_step = [&](auto& sys,
                        typename std::vector<scalar>::const_iterator x,
                        std::vector<scalar>& y_end) -> void {
    const std::vector<scalar> y_start(x, x + static_cast<std::ptrdiff_t>(size));
    std::vector<std::vector<scalar>> rate(6, std::vector<scalar>(size));
    std::vector<scalar> stage(size);
    // k1 as the run took it: see Step::step_adjoint.
    bool k1_loaded = false;
    if constexpr (SolvesForValues<System>) {
      if (k1 != nullptr) {
        ode::derivs(sys, y_start, rate[0], k1->time, k1->row);
        k1_loaded = true;
      }
    }
    if (!k1_loaded) {
      ode::derivs(sys, y_start, rate[0], k1 != nullptr ? k1->time : time);
    }
    ++recorded_rates;
    stages(y_start, rate, h, stage,
           [&](int i, const std::vector<scalar>& at,
               std::vector<scalar>& into) -> void {
             if constexpr (SolvesForValues<System>) {
               ode::derivs(sys, at, into, stage_time(i, time, h),
                           std::as_const(solved[i - 1]));
             } else {
               ode::derivs(sys, at, into, stage_time(i, time, h));
             }
             ++recorded_rates;
           });
    step_end(y_start, rate, h, y_end);
  };

  ode::state_and_parameter_adjoints(active, y, lambda_out, whole_step, lambda_in,
                                    parameter_adjoint);
}

} // namespace ode
} // namespace odelia

#endif
