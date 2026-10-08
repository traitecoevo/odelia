// -*-c++-*-
#ifndef ODELIA_ODE_STEADY_STATE_HPP_
#define ODELIA_ODE_STEADY_STATE_HPP_

// Steady-state solve + implicit-function-theorem parameter sensitivity for an
// autonomous System.
//
// For a System whose right-hand side is f(y; theta), this:
//   1. Solves the equilibrium  f(y*, theta) = 0  by damped Newton's method,
//      reusing the state Jacobian J = df/dy the implicit stepper uses
//      (ode_jacobian.hpp: the system's own ode_jacobian() hook, or forward-mode
//      AD on a rebind_from() twin) and the dense LU (ode_linalg.hpp).
//   2. Computes the parameter sensitivity of the fixed point by the implicit
//      function theorem,   dy*/dtheta = - (df/dy)^-1 (df/dtheta),   reusing the
//      LU factorization of df/dy already formed at the solution and the
//      forward-AD parameter Jacobian df/dtheta (Jacobian::compute_params).
//   3. Reports whether the fixed point is attracting from the eigenvalues of
//      df/dy (all in the left half-plane), computed from the same matrix.
//
// This is deliberately *endpoint-only*: no time integration through the
// transient, and no nested AD. Everything runs at a passive scalar type; the
// Jacobians use one tape-free forward-mode layer internally. The sensitivity
// comes back as plain rows, so a caller whose own model sits on an adjoint tape
// attaches them to y* as a supplied derivative (implicit_node.hpp) rather
// than taping the Newton iteration, which is why solve() refuses an active
// scalar type outright. An optional warm-start integrates the transient to
// reach the attracting basin before Newton.
//
// Newton finds *a* root of f, not necessarily an attracting one: from a guess
// near a saddle or a repeller (the trivial equilibrium of a demographic model,
// say) it converges there in a step or two. solve() reports what it found and
// is_stable() says which kind it is; solve_with_warmup() is the call that seeks
// an attractor, and integrates the transient whenever the Newton root is not
// one.
//
// Scope: fixed-point attractors of autonomous systems only. Equilibrium is not
// well-defined when f depends explicitly on time; time_dependence() surfaces a
// nonzero df/dt at the solution as a guard. A system that declares
// ode_autonomous() is not asked for df/dt at all, as with the implicit stepper.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <type_traits>
#include <vector>
#include <XAD/XAD.hpp>
#include <odelia/ode_util.hpp>
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

  // Requires a state Jacobian from either route (an ode_jacobian() hook, or a
  // rebind_from() hook + non-active scalar). Parameter sensitivity needs the AD
  // route plus the ad_parameters() hook. Callers gate on these; the class
  // instantiates regardless.
  static constexpr bool supported = Jacobian<System>::supported;
  static constexpr bool params_supported = Jacobian<System>::params_supported;

  struct Options {
    double tol = 1e-10;      // convergence: ||f(y)||_inf < tol
    size_t max_iter = 100;   // Newton iteration cap
    bool line_search = true; // backtracking on ||f||_2 for a wider basin
    double min_lambda = 1e-10; // smallest line-search fraction before giving up
  };

  struct Result {
    state_type y;              // equilibrium estimate
    state_type residual;       // f(y) at the estimate
    double residual_norm = 0;  // ||f(y)||_inf
    size_t iterations = 0;     // Newton iterations taken
    bool converged = false;
    // Not converged with iterations < max_iter means Newton stalled: no step
    // fraction down to min_lambda reduced ||f||, the step was not finite, or
    // df/dy was singular at the iterate.
    bool warmed = false;       // solve_with_warmup() integrated the transient
  };

  void resize(size_t n_) {
    n = n_;
    jac.resize(n_);
  }

  // Solve f(y) = 0 by damped Newton starting from `y0`. The evaluation time is
  // the system's current autonomous time (0 for time-homogeneous systems). On
  // convergence the state Jacobian df/dy at the solution is retained and
  // factored for sensitivity() and the eigenvalue queries, which refuse to
  // answer otherwise. The system is left set at the final estimate.
  Result solve(System& system, const state_type& y0, const Options& opt = Options()) {
    static_assert(std::is_arithmetic<value_type>::value,
                  "SteadyState runs at a passive scalar type: taping a Newton "
                  "iteration is never the right derivative. Put y* on an outer "
                  "tape by attaching sensitivity() as a supplied derivative.");
    if constexpr (!supported) {
      util::stop("Steady-state solve needs a Jacobian: the system must provide "
                 "an ode_jacobian() hook, or a rebind_from() hook and a non-active "
                 "scalar type.");
      return Result();
    } else {
      resize(system.ode_size());
      t_eval = ode::ode_time(system);
      solved = false;
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

        // Newton direction: solve (df/dy) dy = -f, factoring df/dy afresh. `f`
        // is f(y) already in hand, which a hook (finite differences) can reuse.
        // A singular df/dy at the iterate is a refusal from the LU; here
        // it means Newton cannot proceed from this point, so the solve stops
        // unconverged and a warm start (below) can take over.
        jac.compute(system, y, t_eval, f, J);
        lu = J;
        try {
          linalg::lu_decompose(lu, n, piv);
        } catch (const util::DomainError&) {
          break;
        }
        for (size_t i = 0; i < n; ++i) {
          neg_f[i] = -f[i];
        }
        linalg::lu_solve(lu, n, piv, neg_f, dy);

        // Backtracking line search on ||f||_2 (Armijo-style) to widen the basin
        // of convergence; a full Newton step is tried first. A step is taken
        // only if it is finite and (with the line search on) reduces ||f||; if
        // no fraction down to min_lambda does, the iteration has stalled and
        // the solve stops rather than taking a step that makes things worse.
        double lambda = 1.0;
        bool accepted = false;
        while (true) {
          for (size_t i = 0; i < n; ++i) {
            y_try[i] = y[i] + lambda * dy[i];
          }
          ode::derivs(system, y_try, f_try, t_eval);
          const double ftnorm2 = norm2(f_try);
          const bool decreased =
              !opt.line_search || ftnorm2 < (1.0 - 1e-4 * lambda) * fnorm2;
          if (std::isfinite(ftnorm2) && decreased) {
            y = y_try;
            f = f_try;
            fnorm2 = ftnorm2;
            accepted = true;
            break;
          }
          if (lambda <= opt.min_lambda) {
            break;
          }
          lambda *= 0.5;
        }
        if (!accepted) {
          break;
        }
      }

      if (norm_inf(f) < opt.tol) {
        res.converged = true;
      }
      res.y = y;
      res.residual = f;
      res.residual_norm = norm_inf(f);
      y_star = y;

      if (res.converged) {
        // Form df/dy once at the solution for the eigenvalue/stability queries
        // and factor it for sensitivity(). A singular df/dy at a genuine root
        // is a bifurcation point: the eigenvalues still answer (one of them is
        // zero) but the implicit function theorem does not apply, so only the
        // factorisation is withheld.
        jac.compute(system, y_star, t_eval, f, J);
        solved = true;
        lu = J;
        try {
          linalg::lu_decompose(lu, n, piv);
          factored = true;
        } catch (const util::DomainError&) {
          factored = false;
        }

        // Record ||df/dt|| as an autonomy guard, unless the system declares
        // itself autonomous.
        if (ode::is_autonomous(system)) {
          dfdt_norm = 0.0;
        } else {
          state_type dfdt;
          dfdt_fd(system, y_star, t_eval, f, dfdt);
          dfdt_norm = norm_inf(dfdt);
        }
      }

      // The Jacobian hook and the df/dt difference both evaluate the system
      // elsewhere; put it back on the estimate before handing it back.
      ode::derivs(system, y_star, f, t_eval);
      return res;
    }
  }

  // Seek an *attracting* fixed point. Newton is tried from `y0` first; if it
  // does not converge, or converges to a root that is not attracting (a saddle
  // or repeller, which Newton finds as readily as an attractor), the transient
  // is integrated with the requested method (RODAS for a stiff system) over
  // `times` to reach the attracting basin, and Newton is retried from the
  // integrated endpoint. `warmed` in the Result says whether that happened.
  // A root reached this way has been approached dynamically, which is as much
  // confirmation of attraction as the spectral abscissa gives, and holds when
  // the abscissa is too close to zero to read.
  Result solve_with_warmup(System& system, const state_type& y0,
                           const OdeControl& control, Method method,
                           const std::vector<double>& times,
                           const Options& opt = Options()) {
    Result res = solve(system, y0, opt);
    if (times.size() < 2 || (res.converged && is_stable())) {
      return res;
    }

    // Integrate a copy of the system forward from y0 to approach equilibrium.
    Solver<System> solver(system, control, method);
    std::vector<double> y0d(y0.size());
    for (size_t i = 0; i < y0.size(); ++i) {
      y0d[i] = util::to_passive(y0[i]);
    }
    solver.set_collect(false);
    solver.set_state(y0d, times.front());
    solver.advance_adaptive(times);

    // The solver stepped a copy. Hand it back: its evaluation counters and
    // whatever it cached on the way are the system Newton continues from.
    system = solver.get_system_ref();
    res = solve(system, solver.state(), opt);
    res.warmed = true;
    return res;
  }

  // Parameter sensitivity of the fixed point via the implicit function theorem,
  //   dy*/dtheta = - (df/dy)^-1 (df/dtheta),
  // returned row-major as an n x n_params matrix, S[row*n_params + col] =
  // d y*_row / d theta_col, reusing the retained LU factorization of df/dy.
  // `n_params` is set to the number of parameters the System exposes. Refuses
  // unless the last solve() converged and df/dy was nonsingular there.
  //
  // The columns are only as complete as the System's ad_parameters() contract
  // (ode_jacobian.hpp): a parameter the rates read through a quantity cached at
  // construction contributes a zero column, silently. check_parameters() below
  // detects that.
  std::vector<value_type> sensitivity(const System& system, size_t& n_params) {
    n_params = 0;
    if (!solved) {
      util::stop("sensitivity() needs a converged solve() first.");
    }
    if (!factored) {
      util::stop("df/dy is singular at the equilibrium (a bifurcation point): "
                 "the implicit function theorem does not give a sensitivity.");
    }
    if constexpr (!params_supported) {
      util::stop("Parameter sensitivity needs an ad_parameters() hook on the "
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

  // Self-check of the ad_parameters() contract at the equilibrium: the largest
  // relative discrepancy, over every entry, between the forward-AD df/dtheta
  // and a central difference of f taken through rebind_from(). The difference
  // perturbs a parameter on a copy of the System and then *rebuilds* the System
  // from that copy's values, so anything the constructor derives from the
  // parameter is derived afresh; the AD column, which seeds the parameter in
  // place, sees none of that. A System whose rates read every parameter live
  // agrees to finite-difference accuracy (~1e-6 at the default step); one with
  // a cached derived quantity returns a discrepancy of order one. Costs two
  // rebinds and two rate evaluations per parameter.
  double check_parameters(const System& system, double rel_step = 1e-6) {
    if (!solved) {
      util::stop("check_parameters() needs a converged solve() first.");
    }
    if constexpr (!params_supported) {
      util::stop("check_parameters() needs an ad_parameters() hook on the "
                 "system exposing the differentiable parameters.");
      return 0.0;
    } else {
      size_t n_params = 0;
      std::vector<value_type> Jp;
      jac.compute_params(system, y_star, t_eval, Jp, n_params);

      System probe = system.template rebind_from<value_type>();
      std::vector<value_type*> params = probe.ad_parameters();
      state_type f_plus(n), f_minus(n);
      double worst = 0.0;
      for (size_t c = 0; c < n_params; ++c) {
        const value_type p0 = *params[c];
        const double h =
            rel_step * std::max(std::fabs(util::to_passive(p0)), 1.0);
        *params[c] = p0 + value_type(h);
        System up = probe.template rebind_from<value_type>();
        ode::derivs(up, y_star, f_plus, t_eval);
        *params[c] = p0 - value_type(h);
        System down = probe.template rebind_from<value_type>();
        ode::derivs(down, y_star, f_minus, t_eval);
        *params[c] = p0;
        for (size_t r = 0; r < n; ++r) {
          const double fd = util::to_passive(f_plus[r] - f_minus[r]) / (2.0 * h);
          const double ad = util::to_passive(Jp[r * n_params + c]);
          worst = std::max(worst, std::fabs(fd - ad) / (1.0 + std::fabs(fd)));
        }
      }
      return worst;
    }
  }

  // Eigenvalues of df/dy at the converged equilibrium (real/imag parts).
  void eigenvalues(std::vector<double>& re, std::vector<double>& im) const {
    if (!solved) {
      util::stop("eigenvalues() needs a converged solve() first.");
    }
    linalg::eigenvalues(J, n, re, im);
  }

  // The spectral abscissa (largest eigenvalue real part) of df/dy. The fixed
  // point is asymptotically stable (attracting) iff this is < 0.
  double spectral_abscissa() const {
    if (!solved) {
      util::stop("spectral_abscissa() needs a converged solve() first.");
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
  // Norms of the passive values: convergence and line-search decisions are
  // never themselves differentiated quantities.
  static double norm_inf(const state_type& v) {
    double m = 0.0;
    for (const auto& x : v) {
      m = std::max(m, std::fabs(util::to_passive(x)));
    }
    return m;
  }
  static double norm2(const state_type& v) {
    double s = 0.0;
    for (const auto& x : v) {
      const double xv = util::to_passive(x);
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
  bool solved = false;         // last solve() converged; J is at y*
  bool factored = false;       // ... and J was nonsingular there; lu is valid
  double dfdt_norm = 0.0;
};

} // namespace ode
} // namespace odelia

#endif
