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
// What it has that the Cash-Karp pair (ode_step.hpp) has not is a free
// continuous extension of order 4 (Hairer's dopri5.f, `contd5`): the state
// anywhere inside an accepted step, from the step's own stages, at no further
// evaluation and with an error of the same order as the step's. That is what
// makes dense output honest -- output at times the stepper did not land on,
// still within the tolerance it was asked for -- where a cubic Hermite on the
// endpoints is one order short (measured on Lorenz: 40-100x the integration
// error). See SolverInternal::interpolate and Solver::advance_collect.
//
// Matches the `Step` interface -- resize/order/step and the two traits -- so
// SolverInternal drives it through the same adaptive loop; `dense()` is the
// extra. Stage evaluations use the plain 4-argument derivs: no per-stage cache
// hook, which is RKCK's business (plant's replay).

#include <cstddef>
#include <vector>
#include <odelia/ode_interface.hpp>

namespace odelia {
namespace ode {

template <class System>
class DopriStep {
public:
  using value_type = typename System::value_type;
  using state_type = std::vector<value_type>;

  void resize(size_t size_) {
    size = size_;
    for (int s = 0; s < n_stages; ++s) {
      k[s].assign(size, value_type(0.0));
    }
    ytmp.assign(size, value_type(0.0));
    y0.assign(size, value_type(0.0));
    h_last = 0.0;
  }

  // Order of the propagated solution (the controller uses 1/order).
  size_t order() const { return 5; }

  void step(System& system, double time, double step_size,
            state_type& y, state_type& yerr,
            const state_type& dydt_in, state_type& dydt_out) {
    const double h = step_size;
    y0 = y;
    h_last = h;
    k[0] = dydt_in;

    for (size_t i = 0; i < size; ++i) {
      ytmp[i] = y[i] + h * a21 * k[0][i];
    }
    ode::derivs(system, ytmp, k[1], time + c2 * h);

    for (size_t i = 0; i < size; ++i) {
      ytmp[i] = y[i] + h * (a31 * k[0][i] + a32 * k[1][i]);
    }
    ode::derivs(system, ytmp, k[2], time + c3 * h);

    for (size_t i = 0; i < size; ++i) {
      ytmp[i] = y[i] + h * (a41 * k[0][i] + a42 * k[1][i] + a43 * k[2][i]);
    }
    ode::derivs(system, ytmp, k[3], time + c4 * h);

    for (size_t i = 0; i < size; ++i) {
      ytmp[i] = y[i] + h * (a51 * k[0][i] + a52 * k[1][i] + a53 * k[2][i] +
                            a54 * k[3][i]);
    }
    ode::derivs(system, ytmp, k[4], time + c5 * h);

    for (size_t i = 0; i < size; ++i) {
      ytmp[i] = y[i] + h * (a61 * k[0][i] + a62 * k[1][i] + a63 * k[2][i] +
                            a64 * k[3][i] + a65 * k[4][i]);
    }
    ode::derivs(system, ytmp, k[5], time + h);

    // The fifth-order solution; its weights are the seventh stage's row.
    for (size_t i = 0; i < size; ++i) {
      y[i] = y0[i] + h * (b1 * k[0][i] + b3 * k[2][i] + b4 * k[3][i] +
                          b5 * k[4][i] + b6 * k[5][i]);
    }
    ode::derivs(system, y, k[6], time + h);
    dydt_out = k[6];

    // Difference between the fifth- and fourth-order solutions.
    for (size_t i = 0; i < size; ++i) {
      yerr[i] = h * (e1 * k[0][i] + e3 * k[2][i] + e4 * k[3][i] +
                     e5 * k[4][i] + e6 * k[5][i] + e7 * k[6][i]);
    }
  }

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
  static const int n_stages = 7;
  size_t size = 0;
  state_type k[n_stages];
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

} // namespace ode
} // namespace odelia

#endif
