/* Steady state of an R-level right-hand side (#39).
 *
 * Newton on ode::CallbackSystem fed R closures through src/r_system.h, with the
 * Jacobian from the system's hook -- jac(t, y, parms) when given, forward
 * differences through rhs otherwise -- and the stability off its eigenvalues.
 * Passive only. There is no AD route through R code, so the parameter
 * sensitivity is not taken here: ode_steady_state() in R/steady-state.R forms
 * df/dparms by central differences through func and applies the implicit
 * function theorem with the Jacobian returned from this export.
 */

#include <Rcpp.h>
#include <odelia/ode_solver.hpp>
#include <odelia/ode_steady_state.hpp>
#include <odelia/rcpp_interface_helpers.hpp>
#include "r_system.h"

using namespace odelia;

// Solve for the equilibrium of rhs(t0, y, parms) = 0 from y0 by damped Newton,
// integrating the transient over `warmup_times` first (when it has at least
// two entries) whenever Newton's root is not attracting. Returns the estimate
// and its diagnostics; the Jacobian, eigenvalues and stability verdict only
// when the solve converged.
// [[Rcpp::export]]
Rcpp::List RSteadyState_solve(Rcpp::Function rhs, Rcpp::Nullable<Rcpp::Function> jac,
                              SEXP parms, Rcpp::NumericVector y0, double t0,
                              bool autonomous, double jac_fd_step, double jac_fd_floor,
                              double tol, int max_iter, bool line_search,
                              double min_lambda, Rcpp::NumericVector warmup_times,
                              SEXP control_xp, std::string method) {
  using SS = ode::SteadyState<ode::CallbackSystem>;

  std::vector<double> y(y0.begin(), y0.end());
  ode::CallbackSystem sys = rinterface::make_r_system(
      rhs, jac, Rcpp::Nullable<Rcpp::Function>(R_NilValue), Rcpp::RObject(parms), y, t0,
      autonomous, jac_fd_step, jac_fd_floor);

  SS ss;
  SS::Options opt;
  opt.tol = tol;
  opt.max_iter = static_cast<size_t>(max_iter);
  opt.line_search = line_search;
  opt.min_lambda = min_lambda;

  SS::Result res;
  if (warmup_times.size() >= 2) {
    Rcpp::XPtr<ode::OdeControl> ctrl(control_xp);
    std::vector<double> times(warmup_times.begin(), warmup_times.end());
    res = ss.solve_with_warmup(sys, y, *ctrl, parse_method(method), times, opt);
  } else {
    res = ss.solve(sys, y, opt);
  }

  const size_t n = y.size();
  Rcpp::RObject jacobian = R_NilValue;
  Rcpp::RObject eigenvalues = R_NilValue;
  Rcpp::RObject abscissa = R_NilValue;
  Rcpp::RObject stable = R_NilValue;
  if (res.converged) {
    const std::vector<double>& J = ss.state_jacobian();
    Rcpp::NumericMatrix Jm(static_cast<int>(n), static_cast<int>(n));
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        Jm(i, j) = J[i * n + j]; // d f_i / d y_j, as jac(t, y) returns it
      }
    }
    std::vector<double> re, im;
    ss.eigenvalues(re, im);
    Rcpp::ComplexVector ev(static_cast<int>(n));
    for (size_t i = 0; i < n; ++i) {
      ev[i].r = re[i];
      ev[i].i = im[i];
    }
    const double sa = ss.spectral_abscissa();
    jacobian = Jm;
    eigenvalues = ev;
    abscissa = Rcpp::wrap(sa);
    stable = Rcpp::wrap(sa < 0.0);
  }

  return Rcpp::List::create(
      Rcpp::Named("y") = res.y,
      Rcpp::Named("residual") = res.residual,
      Rcpp::Named("residual_norm") = res.residual_norm,
      Rcpp::Named("iterations") = static_cast<int>(res.iterations),
      Rcpp::Named("converged") = res.converged,
      Rcpp::Named("warmed") = res.warmed,
      Rcpp::Named("jacobian") = jacobian,
      Rcpp::Named("eigenvalues") = eigenvalues,
      Rcpp::Named("spectral_abscissa") = abscissa,
      Rcpp::Named("stable") = stable,
      Rcpp::Named("time_dependence") = res.converged ? ss.time_dependence() : NA_REAL,
      Rcpp::Named("n_rhs") = static_cast<double>(sys.n_rhs),
      Rcpp::Named("n_jac") = static_cast<double>(sys.n_jac));
}
