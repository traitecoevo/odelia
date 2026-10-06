/* Solver over an R-level right-hand side (#62).
 *
 * The system is ode::CallbackSystem from the R-free core, fed R closures through
 * src/r_system.h. Passive (double) only: there is no AD route through R code,
 * so these exports take no `active` argument. The R6 class OdeSolver in
 * R/ode-solver.R is the user-facing surface; ode_solve() is the deSolve-shaped
 * convenience over it.
 */

#include <Rcpp.h>
#include <odelia/ode_solver.hpp>
#include <odelia/rcpp_interface_helpers.hpp>
#include "r_system.h"

using namespace odelia;

typedef ode::Solver<ode::CallbackSystem> RSolver;

static Rcpp::XPtr<RSolver> get_rsolver(SEXP xp) {
  return Rcpp::XPtr<RSolver>(xp);
}

// Constructing the solver evaluates rhs once, at (t0, y0), to seed itself.
// With a non-NULL `parms` the callbacks are called as f(t, y, parms).
// [[Rcpp::export]]
SEXP RSolver_new(Rcpp::Function rhs, Rcpp::Nullable<Rcpp::Function> jac,
                 Rcpp::Nullable<Rcpp::Function> state_valid, SEXP parms,
                 Rcpp::NumericVector y0, double t0, SEXP control_xp,
                 std::string method, bool autonomous, double jac_fd_step,
                 double jac_fd_floor) {
  Rcpp::XPtr<ode::OdeControl> ctrl(control_xp);
  std::vector<double> y(y0.begin(), y0.end());
  // A non-positive floor means "the control's absolute tolerance".
  if (!(jac_fd_floor > 0.0)) {
    jac_fd_floor = ctrl->get_tol_abs();
  }
  ode::CallbackSystem sys = rinterface::make_r_system(
      rhs, jac, state_valid, Rcpp::RObject(parms), y, t0, autonomous, jac_fd_step,
      jac_fd_floor);
  auto* solver = new RSolver(sys, *ctrl, parse_method(method));
  // No history: a consumer driving this a step at a time reads the state as it
  // goes, and a copy of the system per step would be pure growth.
  solver->set_collect(false);
  return Rcpp::XPtr<RSolver>(solver, true);
}

// [[Rcpp::export]]
void RSolver_step(SEXP solver_xp, double time_max) {
  get_rsolver(solver_xp)->step(time_max);
}

// [[Rcpp::export]]
void RSolver_advance_adaptive(SEXP solver_xp, Rcpp::NumericVector times) {
  std::vector<double> ts(times.begin(), times.end());
  get_rsolver(solver_xp)->advance_adaptive(ts);
}

// Advance to each of `times` in turn and return the state at each as the rows
// of a matrix, time in the first column: the whole output loop of ode_solve()
// in one call. With `dense` the requested times are read off the interpolant
// of the step that spans them (#24); otherwise each is landed on.
// [[Rcpp::export]]
Rcpp::NumericMatrix RSolver_advance_collect(SEXP solver_xp, Rcpp::NumericVector times,
                                            bool dense) {
  auto solver = get_rsolver(solver_xp);
  std::vector<double> ts(times.begin(), times.end());
  const std::vector<std::vector<double>> rows = solver->advance_collect(ts, dense);
  const size_t n = solver->get_system_ref().ode_size();
  Rcpp::NumericMatrix out(static_cast<int>(rows.size()), static_cast<int>(n + 1));
  for (size_t i = 0; i < rows.size(); ++i) {
    out(i, 0) = ts[i];
    for (size_t j = 0; j < n; ++j) {
      out(i, j + 1) = rows[i][j];
    }
  }
  return out;
}

// [[Rcpp::export]]
double RSolver_time(SEXP solver_xp) {
  return get_rsolver(solver_xp)->time();
}

// [[Rcpp::export]]
Rcpp::NumericVector RSolver_state(SEXP solver_xp) {
  return Rcpp::wrap(get_rsolver(solver_xp)->state());
}

// The rates the system holds for its current state: no evaluation after a step
// or a set_state, which leave it evaluated.
// [[Rcpp::export]]
Rcpp::NumericVector RSolver_rates(SEXP solver_xp) {
  return Rcpp::wrap(ode::r_ode_rates(get_rsolver(solver_xp)->get_system_ref()));
}

// [[Rcpp::export]]
Rcpp::NumericVector RSolver_times(SEXP solver_xp) {
  return Rcpp::wrap(get_rsolver(solver_xp)->times());
}

// Re-seed the state, at any length: the system is resized to y first, then the
// solver reset onto it. Resets the step size to the control's initial value and
// the recorded times; costs one rhs evaluation.
// [[Rcpp::export]]
void RSolver_set_state(SEXP solver_xp, Rcpp::NumericVector y, double time) {
  auto solver = get_rsolver(solver_xp);
  std::vector<double> yy(y.begin(), y.end());
  solver->get_system_ref().resize(yy.size());
  solver->set_state(yy, time);
}

// True after an error in a callback escaped from inside a step: the solver
// holds a half-finished attempt until set_state().
// [[Rcpp::export]]
bool RSolver_mid_step(SEXP solver_xp) {
  return get_rsolver(solver_xp)->mid_step();
}

// [[Rcpp::export]]
double RSolver_step_size(SEXP solver_xp) {
  return get_rsolver(solver_xp)->get_step_size();
}

// [[Rcpp::export]]
void RSolver_set_step_size(SEXP solver_xp, double h) {
  get_rsolver(solver_xp)->set_step_size(h);
}

// Return `value` from the R function whose frame environment is `env`, from
// wherever this is called below it: the non-local exit domain_error() uses to
// leave a callback (see src/r_system.h). `return` is evaluated in that frame
// straight through Rf_eval, under Rcpp's unwind-protect, so the jump to the
// frame is caught as a C++ exception here and resumed at the package boundary.
// Never returns normally.
// [[Rcpp::export(rng = false)]]
SEXP odelia_return_from(SEXP env, SEXP value) {
  if (TYPEOF(env) != ENVSXP) {
    Rcpp::stop("env must be an environment");
  }
  Rcpp::Shield<SEXP> call(Rf_lang2(Rf_install("return"), value));
  return Rcpp::Rcpp_fast_eval(call, env);
}

// [[Rcpp::export]]
Rcpp::NumericVector RSolver_counts(SEXP solver_xp) {
  auto solver = get_rsolver(solver_xp);
  const ode::CallbackSystem& sys = solver->get_system_ref();
  const std::vector<double> times = solver->times();
  const double n_steps = times.empty() ? 0.0 : static_cast<double>(times.size() - 1);
  return Rcpp::NumericVector::create(
      Rcpp::Named("n_rhs") = static_cast<double>(sys.n_rhs),
      Rcpp::Named("n_jac") = static_cast<double>(sys.n_jac),
      Rcpp::Named("n_steps") = n_steps,
      Rcpp::Named("n_rejections") = static_cast<double>(solver->get_n_rejections()));
}
