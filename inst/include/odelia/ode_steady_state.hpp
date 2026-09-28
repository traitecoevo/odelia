// -*-c++-*-
#ifndef ODELIA_ODE_STEADY_STATE_HPP_
#define ODELIA_ODE_STEADY_STATE_HPP_

// Steady-state solve + implicit-function-theorem parameter sensitivity for an
// autonomous System (issue #39, "Case B" of #36).
//
// For a System whose right-hand side is f(y; theta), this:
//   1. Solves the equilibrium  f(y*, theta) = 0  by damped Newton's method,
//      reusing the exact forward-AD state Jacobian J = df/dy (ode_jacobian.hpp)
//      and the dense LU (ode_linalg.hpp).
//   2. Computes the parameter sensitivity of the fixed point by the implicit
//      function theorem,   dy*/dtheta = - (df/dy)^-1 (df/dtheta),   reusing the
//      LU factorization of df/dy already formed for the final Newton step and
//      the forward-AD parameter Jacobian df/dtheta (Jacobian::compute_params).
//   3. Reports whether the fixed point is attracting from the eigenvalues of
//      df/dy (all in the left half-plane), computed from the same matrix.
//
// This is deliberately *endpoint-only*: no time integration through the
// transient, and no nested AD. Everything runs at the solver's scalar type
// (double for the intended passive use); the Jacobians use one tape-free
// forward-mode layer internally. An optional RODAS warm-start integrates the
// transient to reach the attracting basin before Newton, and doubles as a
// dynamical confirmation that the fixed point is attracting.
//
// Scope: fixed-point attractors of autonomous systems only. Equilibrium is not
// well-defined when f depends explicitly on time; time_dependence() surfaces a
// nonzero df/dt at the solution as a guard.

#include <vector>
#include <cstddef>
#include <cmath>
#include <XAD/XAD.hpp>
#include <odelia/ode_interface.hpp>
#include <odelia/ode_jacobian.hpp>
#include <odelia/ode_linalg.hpp>
#include <odelia/ode_solver.hpp>

namespace odelia {
namespace ode {

template <typename System>
class SteadyState {
public:
  using value_type = typename System::value_type;
  using state_type = std::vector<value_type>;

  // Requires the forward-AD state Jacobian (a rebind() hook + non-active
  // scalar). Parameter sensitivity additionally needs the ode_parameters()
  // hook. Callers gate on these; the class instantiates regardless.
  static constexpr bool supported = Jacobian<System>::supported;
  static constexpr bool params_supported = Jacobian<System>::params_supported;

  struct Options {
    double tol = 1e-10;      // convergence: ||f(y)||_inf < tol
    size_t max_iter = 100;   // Newton iteration cap
    bool line_search = true; // backtracking on ||f||_2 for a wider basin
    double min_lambda = 1e-10; // smallest line-search fraction before accepting
  };

  struct Result {
    state_type y;              // equilibrium estimate
    state_type residual;       // f(y) at the estimate
    double residual_norm = 0;  // ||f(y)||_inf
    size_t iterations = 0;     // Newton iterations taken
    bool converged = false;
  };

  void resize(size_t n_) {
    n = n_;
    jac.resize(n_);
  }

  // Solve f(y) = 0 by damped Newton starting from `y0`. The evaluation time is
  // the system's current autonomous time (0 for time-homogeneous systems). On
  // success the state Jacobian df/dy at the solution is factored and retained
  // for sensitivity()/eigenvalue queries; the system is left set at y*.
  Result solve(System& system, const state_type& y0, const Options& opt = Options()) {
    if constexpr (!supported) {
      util::stop("Steady-state solve needs a forward-AD Jacobian: the system "
                 "must provide a rebind() hook and a non-active scalar type.");
      return Result();
    } else {
      resize(system.ode_size());
      t_eval = ode::ode_time(system);
      factored = false;

      state_type y = y0;
      state_type f(n), neg_f(n), dy(n), y_try(n), f_try(n);

      ode::derivs(system, y, f, t_eval);
      double fnorm2 = norm2(f);

      Result res;
      for (size_t iter = 0; iter < opt.max_iter; ++iter) {
        if (norm_inf(f) < opt.tol) {
          res.converged = true;
          break;
        }
        res.iterations = iter + 1;

        // Newton direction: solve (df/dy) dy = -f, factoring df/dy afresh.
        jac.compute(system, y, t_eval, J);
        lu = J;
        linalg::lu_decompose(lu, n, piv);
        for (size_t i = 0; i < n; ++i) {
          neg_f[i] = -f[i];
        }
        linalg::lu_solve(lu, n, piv, neg_f, dy);

        // Backtracking line search on ||f||_2 (Armijo-style) to widen the basin
        // of convergence; a full Newton step is tried first.
        double lambda = 1.0;
        while (true) {
          for (size_t i = 0; i < n; ++i) {
            y_try[i] = y[i] + lambda * dy[i];
          }
          ode::derivs(system, y_try, f_try, t_eval);
          const double ftnorm2 = norm2(f_try);
          const bool accept =
              !opt.line_search ||
              ftnorm2 < (1.0 - 1e-4 * lambda) * fnorm2 ||
              lambda <= opt.min_lambda;
          if (accept) {
            y = y_try;
            f = f_try;
            fnorm2 = ftnorm2;
            break;
          }
          lambda *= 0.5;
        }
      }

      // Settle the system on the final estimate and record the residual.
      ode::derivs(system, y, f, t_eval);
      if (norm_inf(f) < opt.tol) {
        res.converged = true;
      }
      res.y = y;
      res.residual = f;
      res.residual_norm = norm_inf(f);
      y_star = y;

      // Factor df/dy once at the solution for reuse by sensitivity() and the
      // eigenvalue/stability queries. Also record ||df/dt|| as an
      // autonomy/attraction guard.
      jac.compute(system, y_star, t_eval, J);
      lu = J;
      linalg::lu_decompose(lu, n, piv);
      factored = true;

      state_type dfdt;
      dfdt_fd(system, y_star, t_eval, dfdt);
      dfdt_norm = norm_inf(dfdt);

      return res;
    }
  }

  // Try Newton from `y0`; if it fails to converge, integrate the transient with
  // the requested method (RODAS recommended for stiff systems) over `times` to
  // reach the attracting basin, then retry Newton from the integrated endpoint.
  // A successful warm-started solve is also a dynamical confirmation that the
  // fixed point is attracting.
  Result solve_with_warmup(System& system, const state_type& y0,
                           const OdeControl& control, Method method,
                           const std::vector<double>& times,
                           const Options& opt = Options()) {
    Result res = solve(system, y0, opt);
    if (res.converged || times.size() < 2) {
      return res;
    }

    // Integrate a copy of the system forward from y0 to approach equilibrium.
    Solver<System> solver(system, control, method);
    std::vector<double> y0d(y0.size());
    for (size_t i = 0; i < y0.size(); ++i) {
      y0d[i] = xad::value(y0[i]);
    }
    solver.set_collect(false);
    solver.set_state(y0d, times.front());
    solver.advance_adaptive(times);

    return solve(system, solver.state(), opt);
  }

  // Parameter sensitivity of the fixed point via the implicit function theorem,
  //   dy*/dtheta = - (df/dy)^-1 (df/dtheta),
  // returned row-major as an n x n_params matrix, S[row*n_params + col] =
  // d y*_row / d theta_col, reusing the retained LU factorization of df/dy.
  // `n_params` is set to the number of parameters the System exposes.
  std::vector<value_type> sensitivity(System& system, size_t& n_params) {
    n_params = 0;
    if (!factored) {
      util::stop("Call solve() (to convergence) before sensitivity().");
    }
    if constexpr (!params_supported) {
      util::stop("Parameter sensitivity needs an ode_parameters() hook on the "
                 "system exposing the differentiable parameters.");
      return std::vector<value_type>();
    } else {
      std::vector<value_type> Jp; // n x n_params, df/dtheta
      jac.compute_params(system, y_star, t_eval, Jp, n_params);

      std::vector<value_type> S(n * n_params);
      state_type col(n), x(n);
      for (size_t c = 0; c < n_params; ++c) {
        for (size_t r = 0; r < n; ++r) {
          col[r] = -Jp[r * n_params + c]; // right-hand side -df/dtheta_c
        }
        linalg::lu_solve(lu, n, piv, col, x);
        for (size_t r = 0; r < n; ++r) {
          S[r * n_params + c] = x[r];
        }
      }
      return S;
    }
  }

  // Eigenvalues of df/dy at the converged equilibrium (real/imag parts).
  void eigenvalues(std::vector<double>& re, std::vector<double>& im) const {
    if (!factored) {
      util::stop("Call solve() before eigenvalues().");
    }
    linalg::eigenvalues(J, n, re, im);
  }

  // The spectral abscissa (largest eigenvalue real part) of df/dy. The fixed
  // point is asymptotically stable (attracting) iff this is < 0.
  double spectral_abscissa() const {
    if (!factored) {
      util::stop("Call solve() before spectral_abscissa().");
    }
    return linalg::spectral_abscissa(J, n);
  }

  bool is_stable() const { return spectral_abscissa() < 0.0; }

  // ||df/dt||_inf at the solution. Nonzero (beyond finite-difference noise)
  // means f depends explicitly on time, so the "equilibrium" is not a genuine
  // fixed point of an autonomous system -- a guard for misuse on non-autonomous
  // systems, for which equilibrium is not well-defined.
  double time_dependence() const { return dfdt_norm; }

  const state_type& equilibrium() const { return y_star; }
  const std::vector<value_type>& state_jacobian() const { return J; }

private:
  static double norm_inf(const state_type& v) {
    double m = 0.0;
    for (const auto& x : v) {
      m = std::max(m, std::fabs(xad::value(x)));
    }
    return m;
  }
  static double norm2(const state_type& v) {
    double s = 0.0;
    for (const auto& x : v) {
      const double xv = xad::value(x);
      s += xv * xv;
    }
    return std::sqrt(s);
  }

  size_t n = 0;
  double t_eval = 0.0;
  Jacobian<System> jac;
  std::vector<value_type> J;   // df/dy at the solution (row-major n*n)
  std::vector<value_type> lu;  // LU factorization of J
  std::vector<size_t> piv;
  state_type y_star;           // equilibrium estimate
  bool factored = false;
  double dfdt_norm = 0.0;
};

} // namespace ode
} // namespace odelia

#endif
