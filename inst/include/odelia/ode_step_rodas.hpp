// -*-c++-*-
#ifndef ODELIA_ODE_STEP_RODAS_HPP_
#define ODELIA_ODE_STEP_RODAS_HPP_

// RODAS4(3): a 6-stage, L-stable, stiffly-accurate Rosenbrock method of order 4
// with an embedded 3rd-order error estimate (Hairer & Wanner, Solving ODEs II,
// section IV.7; reference implementation `rodas.f`, coefficient set METH=1).
//
// A Rosenbrock (linearly-implicit) method: each stage is one RHS evaluation plus
// one linear solve against a single matrix W = (1/(h*gamma)) I - J, factorized
// once per step. No Newton iteration, so no convergence machinery -- the only
// adaptivity is the existing accuracy-based step-size controller in OdeControl,
// which consumes the embedded error estimate exactly as it does for the explicit
// RKCK stepper.
//
// This class matches the `Step` (RKCK) interface -- resize/order/step and the
// can_use_dydt_in / first_same_as_last traits -- so SolverInternal can drive it
// through the same adaptive loop.
//
// The Jacobian J = df/dy is computed exactly by forward-mode AD (ode_jacobian.hpp),
// which requires the System to expose `template<class U> rebind<U> rebind_from()`. The
// time derivative df/dt (only needed for non-autonomous systems) is a finite
// difference, because the System stores time as a plain double.

#include <vector>
#include <cstddef>
#include <odelia/ode_interface.hpp>
#include <odelia/ode_linalg.hpp>
#include <odelia/ode_jacobian.hpp>

namespace odelia {
namespace ode {

template <class System>
class RodasStep {
public:
  using value_type = typename System::value_type;
  using state_type = std::vector<value_type>;

  // True when the exact-AD Jacobian is instantiable for this scalar type (the
  // passive double solver). False for nested AD types until that path lands.
  static constexpr bool supported = Jacobian<System>::supported;

  void resize(size_t size_) {
    size = size_;
    jac.resize(size_);
    J.assign(size * size, value_type(0.0));
    W.assign(size * size, value_type(0.0));
    dT.assign(size, value_type(0.0));
    arg.assign(size, value_type(0.0));
    ftmp.assign(size, value_type(0.0));
    rhs.assign(size, value_type(0.0));
    for (int s = 0; s < n_stages; ++s) {
      k[s].assign(size, value_type(0.0));
    }
  }

  // Classical order of the method (the controller uses 1/order). The embedded
  // pair is 4(3); we report the higher order 4, matching rodas.f's step exponent.
  size_t order() const { return 4; }

  void step(System& system,
            double time, double step_size,
            state_type& y,
            state_type& yerr,
            const state_type& dydt_in,
            state_type& dydt_out) {
    const double h = step_size;

    // J = df/dy and dT = df/dt, both at the step start (t_n, y_n). W is factored
    // once and reused across all six stage solves.
    jac.compute(system, y, time, J);
    dfdt_fd(system, y, time, dT);

    const value_type fac = value_type(1.0 / (h * gamma));
    for (size_t i = 0; i < size * size; ++i) {
      W[i] = -J[i];
    }
    for (size_t d = 0; d < size; ++d) {
      W[d * size + d] += fac;
    }
    linalg::lu_decompose(W, size, piv);

    // Stage 1: f(t_n, y_n) is already available as dydt_in.
    for (size_t i = 0; i < size; ++i) {
      rhs[i] = dydt_in[i] + h * d1 * dT[i];
    }
    linalg::lu_solve(W, size, piv, rhs, k[0]);

    // Stage 2
    for (size_t i = 0; i < size; ++i) {
      arg[i] = y[i] + a21 * k[0][i];
    }
    ode::derivs(system, arg, ftmp, time + c2 * h);
    for (size_t i = 0; i < size; ++i) {
      rhs[i] = ftmp[i] + (c21 * k[0][i]) / h + h * d2 * dT[i];
    }
    linalg::lu_solve(W, size, piv, rhs, k[1]);

    // Stage 3
    for (size_t i = 0; i < size; ++i) {
      arg[i] = y[i] + a31 * k[0][i] + a32 * k[1][i];
    }
    ode::derivs(system, arg, ftmp, time + c3 * h);
    for (size_t i = 0; i < size; ++i) {
      rhs[i] = ftmp[i] + (c31 * k[0][i] + c32 * k[1][i]) / h + h * d3 * dT[i];
    }
    linalg::lu_solve(W, size, piv, rhs, k[2]);

    // Stage 4
    for (size_t i = 0; i < size; ++i) {
      arg[i] = y[i] + a41 * k[0][i] + a42 * k[1][i] + a43 * k[2][i];
    }
    ode::derivs(system, arg, ftmp, time + c4 * h);
    for (size_t i = 0; i < size; ++i) {
      rhs[i] = ftmp[i] +
               (c41 * k[0][i] + c42 * k[1][i] + c43 * k[2][i]) / h +
               h * d4 * dT[i];
    }
    linalg::lu_solve(W, size, piv, rhs, k[3]);

    // Stage 5 (time node alpha_5 = 1). The stage argument here is the embedded
    // (order-3) solution point.
    for (size_t i = 0; i < size; ++i) {
      arg[i] = y[i] + a51 * k[0][i] + a52 * k[1][i] + a53 * k[2][i] +
               a54 * k[3][i];
    }
    ode::derivs(system, arg, ftmp, time + h);
    for (size_t i = 0; i < size; ++i) {
      rhs[i] = ftmp[i] +
               (c51 * k[0][i] + c52 * k[1][i] + c53 * k[2][i] +
                c54 * k[3][i]) / h;
    }
    linalg::lu_solve(W, size, piv, rhs, k[4]);

    // Stage 6 (time node alpha_6 = 1). arg becomes the embedded solution
    // (previous arg + k5); the full order-4 solution is arg + k6.
    for (size_t i = 0; i < size; ++i) {
      arg[i] += k[4][i];
    }
    ode::derivs(system, arg, ftmp, time + h);
    for (size_t i = 0; i < size; ++i) {
      rhs[i] = ftmp[i] +
               (c61 * k[0][i] + c62 * k[1][i] + c63 * k[2][i] +
                c64 * k[3][i] + c65 * k[4][i]) / h;
    }
    linalg::lu_solve(W, size, piv, rhs, k[5]);

    // Update and embedded error estimate. y_{n+1} = arg + k6; the difference
    // between the order-4 and order-3 solutions is exactly k6.
    for (size_t i = 0; i < size; ++i) {
      y[i] = arg[i] + k[5][i];
      yerr[i] = k[5][i];
    }

    // Derivative at the new point, for the controller and FSAL bookkeeping.
    ode::derivs(system, y, dydt_out, time + h);
  }

  // A single J and factorization are used for the whole step; the start-of-step
  // derivative (dydt_in) is genuinely f(t_n, y_n) and is reused for stage 1.
  static const bool can_use_dydt_in = true;
  // RODAS is stiffly accurate but we recompute dydt_out explicitly, so we do not
  // claim FSAL reuse of dydt_out as the next dydt_in.
  static const bool first_same_as_last = false;

private:
  static const int n_stages = 6;

  size_t size = 0;
  Jacobian<System> jac;
  std::vector<value_type> J;   // row-major n*n, df/dy
  std::vector<value_type> W;   // row-major n*n, (1/(h*gamma)) I - J
  std::vector<size_t> piv;
  std::vector<value_type> dT;  // df/dt
  std::vector<value_type> arg; // stage argument / running solution
  std::vector<value_type> ftmp;
  std::vector<value_type> rhs;
  state_type k[n_stages];

  // RODAS coefficients (Hairer & Wanner, rodas.f METH=1).
  static constexpr double gamma = 0.25;
  static constexpr double c2 = 0.386;
  static constexpr double c3 = 0.21;
  static constexpr double c4 = 0.63;
  static constexpr double d1 = 0.25;
  static constexpr double d2 = -0.1043;
  static constexpr double d3 = 0.1035;
  static constexpr double d4 = -0.03620000000000023;
  static constexpr double a21 = 1.544;
  static constexpr double a31 = 0.9466785280815826;
  static constexpr double a32 = 0.2557011698983284;
  static constexpr double a41 = 3.314825187068521;
  static constexpr double a42 = 2.896124015972201;
  static constexpr double a43 = 0.9986419139977817;
  static constexpr double a51 = 1.221224509226641;
  static constexpr double a52 = 6.019134481288629;
  static constexpr double a53 = 12.53708332932087;
  static constexpr double a54 = -0.6878860361058950;
  static constexpr double c21 = -5.6688;
  static constexpr double c31 = -2.430093356833875;
  static constexpr double c32 = -0.2063599157091915;
  static constexpr double c41 = -0.1073529058151375;
  static constexpr double c42 = -9.594562251023355;
  static constexpr double c43 = -20.47028614809616;
  static constexpr double c51 = 7.496443313967647;
  static constexpr double c52 = -10.24680431464352;
  static constexpr double c53 = -33.99990352819905;
  static constexpr double c54 = 11.70890893206160;
  static constexpr double c61 = 8.083246795921522;
  static constexpr double c62 = -7.981132988064893;
  static constexpr double c63 = -31.52159432874371;
  static constexpr double c64 = 16.31930543123136;
  static constexpr double c65 = -6.058818238834054;
};

} // namespace ode
} // namespace odelia

#endif
