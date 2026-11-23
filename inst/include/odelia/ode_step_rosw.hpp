// -*-c++-*-
#ifndef ODELIA_ODE_STEP_ROSW_HPP_
#define ODELIA_ODE_STEP_ROSW_HPP_

// ROS34PW2: a 4-stage, L-stable, stiffly-accurate Rosenbrock-W method of order
// 3 with an embedded 2nd-order error estimate (J. Rang and L. Angermann, New
// Rosenbrock W-methods of order 3 for partial differential algebraic equations
// of index 1, BIT Numer. Math. 45 (2005) 761-787). The tableau below is the
// paper's.
//
// A W-method's order conditions hold for any matrix in place of the Jacobian,
// and so do its embedded pair's, so the error estimate stays honest when the
// matrix is out of date. That is what lets this stepper keep J across steps:
// it is formed at a step's start and then serves the steps after it, until a
// step it served is rejected or it reaches the age the control allows
// (OdeControl::jacobian_max_age). W = (1/(h*gamma)) I - J is formed and
// factored at every step, at that step's h, so a change of step size never
// makes the kept J wrong. Where the Jacobian is cheap, set the age to 1: J is
// then formed at every step and the method is an ordinary Rosenbrock method of
// order 3.
//
// df/dt is formed and kept with J: the time row of the autonomous form's
// Jacobian is part of the matrix a W-method may approximate.
//
// The stages are solved against W in the form of Hairer and Wanner, Solving
// ODEs II, section IV.7, eq. 7.25, as RODAS's are, with the coefficients that
// form takes derived below from the paper's alpha, Gamma, b and b-hat.

#include <array>
#include <vector>
#include <cstddef>
#include <odelia/ode_interface.hpp>
#include <odelia/ode_linalg.hpp>
#include <odelia/ode_jacobian.hpp>

namespace odelia {
namespace ode {

namespace rosw {

constexpr std::size_t s = 4;
using matrix = std::array<std::array<double, s>, s>;
using row = std::array<double, s>;

constexpr double gamma = 4.3586652150845900e-01;
constexpr matrix alpha{{
  {0, 0, 0, 0},
  {8.7173304301691801e-01, 0, 0, 0},
  {8.4457060015369423e-01, -1.1299064236484185e-01, 0, 0},
  {0, 0, 1, 0}}};
constexpr matrix Gamma{{
  {gamma, 0, 0, 0},
  {-8.7173304301691801e-01, gamma, 0, 0},
  {-9.0338057013044082e-01, 5.4180672388095326e-02, gamma, 0},
  {2.4212380706095346e-01, -1.2232505839045147e+00, 5.4526025533510214e-01,
   gamma}}};
constexpr row b{2.4212380706095346e-01, -1.2232505839045147e+00,
                1.5452602553351020e+00, 4.3586652150845900e-01};
constexpr row b_hat{3.7810903145819369e-01, -9.6042292212423178e-02,
                    5.0000000000000000e-01, 2.1793326075422950e-01};

// Gamma^{-1}, lower triangular, by forward substitution.
constexpr matrix inverse(const matrix& g) {
  matrix r{};
  for (std::size_t j = 0; j < s; ++j) {
    r[j][j] = 1.0 / g[j][j];
    for (std::size_t i = j + 1; i < s; ++i) {
      double acc = 0.0;
      for (std::size_t k = j; k < i; ++k) {
        acc += g[i][k] * r[k][j];
      }
      r[i][j] = -acc / g[i][i];
    }
  }
  return r;
}
constexpr matrix Gamma_inv = inverse(Gamma);

constexpr matrix times(const matrix& x, const matrix& y) {
  matrix r{};
  for (std::size_t i = 0; i < s; ++i) {
    for (std::size_t j = 0; j < s; ++j) {
      for (std::size_t k = 0; k < s; ++k) {
        r[i][j] += x[i][k] * y[k][j];
      }
    }
  }
  return r;
}
constexpr row times(const row& x, const matrix& y) {
  row r{};
  for (std::size_t j = 0; j < s; ++j) {
    for (std::size_t k = 0; k < s; ++k) {
      r[j] += x[k] * y[k][j];
    }
  }
  return r;
}

// The transformed coefficients of eq. 7.25: a = alpha Gamma^{-1} for the stage
// arguments, c = diag(1/gamma) - Gamma^{-1} for the stage right-hand sides, m
// and m_hat = b and b_hat times Gamma^{-1} for the solution and the embedded
// one; and each stage's time node and df/dt coefficient, the row sums of alpha
// and Gamma.
constexpr matrix a = times(alpha, Gamma_inv);
constexpr matrix c = [] {
  matrix r{};
  for (std::size_t i = 0; i < s; ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      r[i][j] = -Gamma_inv[i][j];
    }
  }
  return r;
}();
constexpr row m = times(b, Gamma_inv);
constexpr row m_hat = times(b_hat, Gamma_inv);
constexpr row row_sum(const matrix& x) {
  row r{};
  for (std::size_t i = 0; i < s; ++i) {
    for (std::size_t j = 0; j < s; ++j) {
      r[i] += x[i][j];
    }
  }
  return r;
}
constexpr row node = row_sum(alpha);
constexpr row d = row_sum(Gamma);

} // namespace rosw

template <class System>
class RosWStep {
public:
  using value_type = typename System::value_type;
  using state_type = std::vector<value_type>;

  static constexpr bool supported = Jacobian<System>::supported;

  // J and df/dt, where they were formed, how many steps they have served, and
  // where the last attempt started. All of it decides the next step's
  // Jacobian, so a caller taking a step back puts all of it back.
  struct kept_jacobian {
    std::vector<value_type> J;  // row-major n*n, df/dy
    std::vector<value_type> dT; // df/dt
    bool valid = false;
    double t = 0.0;
    std::vector<value_type> y;
    std::size_t uses = 0;
    bool attempted = false;
    double t_attempt = 0.0;
    std::vector<value_type> y_attempt;
  };

  void resize(size_t size_) {
    size = size_;
    jac.resize(size_);
    // A kept Jacobian belongs to a state of this size and to the run that
    // formed it; the solver resizes on every load of a state, so a new run
    // starts without one.
    kept = kept_jacobian{};
    kept.dT.assign(size, value_type(0.0));
    arg.assign(size, value_type(0.0));
    ftmp.assign(size, value_type(0.0));
    rhs.assign(size, value_type(0.0));
    for (auto& u_i : u) {
      u_i.assign(size, value_type(0.0));
    }
  }

  // Order of the solution; the embedded estimate is order 2.
  size_t order() const { return 3; }

  // `max_age`, the most steps one J may serve; 0 is no limit.
  void step(System& system,
            double time, double step_size,
            state_type& y,
            state_type& yerr,
            const state_type& dydt_in,
            state_type& dydt_out,
            std::size_t max_age) {
    const double h = step_size;
    using rosw::s;

    // An attempt from where the last one started is a retry of a rejected
    // step (the solver restores y by copy, so the comparison is exact). A
    // retry keeps J only where J was formed at this point: a J from an
    // earlier point may be why the step was rejected.
    const bool retry = kept.attempted && time == kept.t_attempt &&
                       y == kept.y_attempt;
    bool form;
    if (retry) {
      form = !(kept.valid && time == kept.t && y == kept.y);
    } else {
      form = !kept.valid || (max_age > 0 && kept.uses >= max_age);
      if (!form) {
        ++kept.uses;
      }
    }
    if (form) {
      jac.compute(system, y, time, dydt_in, kept.J);
      if (ode::is_autonomous(system)) {
        kept.dT.assign(size, value_type(0.0));
      } else {
        dfdt_fd(system, y, time, dydt_in, kept.dT);
      }
      kept.valid = true;
      kept.t = time;
      kept.y = y;
      kept.uses = 1;
    }
    kept.attempted = true;
    kept.t_attempt = time;
    kept.y_attempt = y;

    const std::vector<value_type>& J = kept.J;
    const std::vector<value_type>& dT = kept.dT;
    W.resize(size * size);
    const value_type fac = value_type(1.0 / (h * rosw::gamma));
    for (size_t i = 0; i < size * size; ++i) {
      W[i] = -J[i];
    }
    for (size_t q = 0; q < size; ++q) {
      W[q * size + q] += fac;
    }
    linalg::lu_decompose(W, size, piv);

    // Stage 1 is at the step's start, where dydt_in is f.
    for (size_t q = 0; q < size; ++q) {
      rhs[q] = dydt_in[q] + h * rosw::d[0] * dT[q];
    }
    linalg::lu_solve(W, size, piv, rhs, u[0]);

    for (size_t i = 1; i < s; ++i) {
      for (size_t q = 0; q < size; ++q) {
        value_type acc = y[q];
        for (size_t j = 0; j < i; ++j) {
          acc += rosw::a[i][j] * u[j][q];
        }
        arg[q] = acc;
      }
      ode::derivs(system, arg, ftmp, time + rosw::node[i] * h);
      for (size_t q = 0; q < size; ++q) {
        value_type acc = 0.0;
        for (size_t j = 0; j < i; ++j) {
          acc += rosw::c[i][j] * u[j][q];
        }
        rhs[q] = ftmp[q] + acc / h + h * rosw::d[i] * dT[q];
      }
      linalg::lu_solve(W, size, piv, rhs, u[i]);
    }

    for (size_t q = 0; q < size; ++q) {
      value_type sol = y[q], err = 0.0;
      for (size_t i = 0; i < s; ++i) {
        sol += rosw::m[i] * u[i][q];
        err += (rosw::m[i] - rosw::m_hat[i]) * u[i][q];
      }
      y[q] = sol;
      yerr[q] = err;
    }

    // Derivative at the new point, for the controller and FSAL bookkeeping.
    ode::derivs(system, y, dydt_out, time + h);
  }

  // What decides the next step's Jacobian, for a caller that takes a step back
  // and must leave the stepper as it found it.
  const kept_jacobian& jacobian_state() const { return kept; }
  void set_jacobian_state(const kept_jacobian& state) { kept = state; }

  // As RODAS: stage 1 uses f(t_n, y_n), and dydt_out, evaluated at the new
  // point, is the next step's f(t_n, y_n).
  static const bool can_use_dydt_in = true;
  static const bool first_same_as_last = true;

private:
  size_t size = 0;
  Jacobian<System> jac;
  kept_jacobian kept;
  std::vector<value_type> W;   // row-major n*n, (1/(h*gamma)) I - J
  std::vector<size_t> piv;
  std::vector<value_type> arg;
  std::vector<value_type> ftmp;
  std::vector<value_type> rhs;
  state_type u[rosw::s];
};

} // namespace ode
} // namespace odelia

#endif
