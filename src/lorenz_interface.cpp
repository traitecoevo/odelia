/* Lorenz system interface for odelia package
 * 
 * This file provides Lorenz-specific functions:
 * - LorenzSystem creation and manipulation
 * - Solver creation (calls generic templates with LorenzSystem type)
 * 
 * Most Solver_* functions are thin wrappers that call generic implementations
 * from ode_interface.cpp with LorenzSystem types.
 */

#include <Rcpp.h>
#include <XAD/XAD.hpp>
#include <odelia/ode_solver.hpp>
#include <odelia/ode_fit.hpp>
#include <odelia/solver_interface.hpp>
#include <examples/lorenz_system.hpp>

using namespace Rcpp;
using namespace odelia;

// Define types for Lorenz system
typedef LorenzSystem<double> SystemType;
typedef LorenzSystem<xad::adj<double>::active_type> ActiveSystemType;

// Helper to get system pointer
inline Rcpp::XPtr<SystemType> get_system(SEXP xp) {
  return Rcpp::XPtr<SystemType>(xp);
}

//-------------------------------------------------------------------------
// System interface (Lorenz-specific)

// [[Rcpp::export]]
SEXP System_new(double sigma, double R, double b) {
  return Rcpp::XPtr<SystemType>(new SystemType(sigma, R, b), true);
}

// [[Rcpp::export]]
Rcpp::NumericVector System_pars(SEXP system_xp) {
  return Rcpp::wrap(get_system(system_xp)->pars());
}

// [[Rcpp::export]]
void System_set_params(SEXP system_xp, Rcpp::NumericVector params) {
  auto lor = get_system(system_xp);
  std::vector<double> tmp(params.begin(), params.end());
  lor->set_params(tmp.begin());
  lor->compute_rates();
}

// [[Rcpp::export]]
void System_set_state(SEXP system_xp, Rcpp::NumericVector y, double time) {
  auto lor = get_system(system_xp);
  if (y.size() != lor->ode_size()) {
    Rcpp::stop("State vector size mismatch");
  }
  std::vector<double> tmp(y.begin(), y.end());
  lor->set_ode_state(tmp.begin(), time);
}

// [[Rcpp::export]]
Rcpp::NumericVector System_state(SEXP system_xp) {
  auto lor = get_system(system_xp);
  std::vector<double> tmp(lor->ode_size());
  lor->ode_state(tmp.begin());
  return Rcpp::wrap(tmp);
}

// [[Rcpp::export]]
void System_set_initial_state(SEXP system_xp, Rcpp::NumericVector y, double t0 = 0.0) {
  auto lor = get_system(system_xp);
  if (y.size() != lor->ode_size()) {
    Rcpp::stop("State vector size mismatch");
  }
  std::vector<double> tmp(y.begin(), y.end());
  lor->set_initial_state(tmp.begin(), t0);
}

// [[Rcpp::export]]
void System_reset(SEXP system_xp) {
  auto lor = get_system(system_xp);
  lor->reset();
}

// [[Rcpp::export]]
Rcpp::NumericVector System_rates(SEXP system_xp) {
  auto lor = get_system(system_xp);
  std::vector<double> tmp(lor->ode_size());
  lor->ode_rates(tmp.begin());
  return Rcpp::wrap(tmp);
}

//-------------------------------------------------------------------------
// Solver interface (Lorenz-specific creation, generic operations)

// Solver creation - Lorenz-specific (must know LorenzSystem type)
// [[Rcpp::export]]
SEXP Solver_new(SEXP system_xp, SEXP control_xp, bool active = false) {
  Rcpp::XPtr<SystemType> sys(system_xp);
  Rcpp::XPtr<ode::OdeControl> ctrl(control_xp);
  
  if (active) {
    // Initialize new AD compatible system
    auto params = sys->pars();
    
    std::vector<double> initial_state(sys->ode_size());
    sys->ode_initial_state(initial_state.begin());
    auto t0 = sys->ode_t0();
    
    ActiveSystemType sys_active(params[0], params[1], params[2]);
    sys_active.set_initial_state(initial_state.begin(), t0);
    
    auto* solver = new ode::Solver<ActiveSystemType>(sys_active, *ctrl);
    
    return Rcpp::XPtr<ode::Solver<ActiveSystemType>>(solver, true);
  } else {
    // Create regular solver
    SystemType sys_copy(*sys);
    auto* solver = new ode::Solver<SystemType>(sys_copy, *ctrl);
    
    return Rcpp::XPtr<ode::Solver<SystemType>>(solver, true);
  }
}

// All other Solver functions call generic templates with LorenzSystem types

// [[Rcpp::export]]
void Solver_reset(SEXP solver_xp, bool active = false) {
  odelia::solver::Solver_reset_impl<SystemType, ActiveSystemType>(solver_xp, active);
}

// [[Rcpp::export]]
double Solver_time(SEXP solver_xp, bool active = false) {
  return odelia::solver::Solver_time_impl<SystemType, ActiveSystemType>(solver_xp, active);
}

// [[Rcpp::export]]
Rcpp::NumericVector Solver_state(SEXP solver_xp, bool active = false) {
  return odelia::solver::Solver_state_impl<SystemType, ActiveSystemType>(solver_xp, active);
}

// [[Rcpp::export]]
Rcpp::NumericVector Solver_times(SEXP solver_xp, bool active = false) {
  return odelia::solver::Solver_times_impl<SystemType, ActiveSystemType>(solver_xp, active);
}

// [[Rcpp::export]]
void Solver_set_state(SEXP solver_xp, Rcpp::NumericVector y, double time, bool active = false) {
  odelia::solver::Solver_set_state_impl<SystemType, ActiveSystemType>(solver_xp, y, time, active);
}

// [[Rcpp::export]]
void Solver_advance_adaptive(SEXP solver_xp, Rcpp::NumericVector times, bool active = false) {
  odelia::solver::Solver_advance_adaptive_impl<SystemType, ActiveSystemType>(solver_xp, times, active);
}

// [[Rcpp::export]]
void Solver_advance_fixed(SEXP solver_xp, Rcpp::NumericVector times, bool active = false) {
  odelia::solver::Solver_advance_fixed_impl<SystemType, ActiveSystemType>(solver_xp, times, active);
}

// [[Rcpp::export]]
void Solver_step(SEXP solver_xp, bool active = false) {
  odelia::solver::Solver_step_impl<SystemType, ActiveSystemType>(solver_xp, active);
}

// [[Rcpp::export]]
bool Solver_get_collect(SEXP solver_xp, bool active = false) {
  return odelia::solver::Solver_get_collect_impl<SystemType, ActiveSystemType>(solver_xp, active);
}

// [[Rcpp::export]]
void Solver_set_collect(SEXP solver_xp, bool x, bool active = false) {
  odelia::solver::Solver_set_collect_impl<SystemType, ActiveSystemType>(solver_xp, x, active);
}

// [[Rcpp::export]]
std::size_t Solver_get_history_size(SEXP solver_xp, bool active = false) {
  return odelia::solver::Solver_get_history_size_impl<SystemType, ActiveSystemType>(solver_xp, active);
}

// Helper to get column names (Lorenz-specific)
CharacterVector get_column_names() {
  return CharacterVector::create("time", "x", "y", "z", "dxdt", "dydt", "dzdt");
}

// [[Rcpp::export]]
Rcpp::DataFrame Solver_get_history_step(SEXP solver_xp, std::size_t i, bool active = false) {
  return odelia::solver::Solver_get_history_step_impl<SystemType, ActiveSystemType>(
    solver_xp, i, get_column_names(), active
  );
}

// [[Rcpp::export]]
Rcpp::List Solver_get_history(SEXP solver_xp, bool active = false) {
  return odelia::solver::Solver_get_history_impl<SystemType, ActiveSystemType>(
    solver_xp, get_column_names(), active
  );
}

// [[Rcpp::export]]
void Solver_set_target(SEXP solver_xp, 
                      Rcpp::NumericVector times,
                      Rcpp::NumericMatrix target,
                      Rcpp::IntegerVector obs_indices,
                      bool active = false) {
  odelia::solver::Solver_set_target_impl<SystemType, ActiveSystemType>(
    solver_xp, times, target, obs_indices, active
  );
}

// [[Rcpp::export]]
Rcpp::List Solver_fit(SEXP solver_xp,
                     Rcpp::Nullable<Rcpp::NumericVector> ic = R_NilValue,
                     Rcpp::Nullable<Rcpp::NumericVector> params = R_NilValue) {
  return odelia::solver::Solver_fit_impl<SystemType, ActiveSystemType>(
    solver_xp, ic, params
  );
}

//-------------------------------------------------------------------------
// Comparison function for deSolve

// [[Rcpp::export]]
List lorenz_rhs(double t, NumericVector state, NumericVector pars) {
  double x = state["x"], y = state["y"], z = state["z"];
  double sigma = pars["sigma"], rho = pars["rho"], beta = pars["beta"];
  
  double dx = sigma * (y - x);
  double dy = x * (rho - z) - y;
  double dz = x * y - beta * z;
  
  return List::create(NumericVector::create(dx, dy, dz));
}

//-------------------------------------------------------------------------
// Test function for parameter type conversion

// [[Rcpp::export]]
Rcpp::List test_param_types(SEXP system_xp) {
  using ad = xad::adj<double>;
  using ad_type = ad::active_type;
  using ActiveSys = LorenzSystem<ad_type>;
  
  Rcpp::XPtr<SystemType> sys(system_xp);
  
  // Get params from passive system (doubles)
  auto params = sys->pars();
  
  // Create active system using double params (should cast to ad_type)
  auto* sys_active = new ActiveSys(params[0], params[1], params[2]);
  
  // Get params back from active system
  auto params_active = sys_active->pars();
  
  // Test if gradient flows through params
  ad::tape_type tape;
  
  std::vector<ad_type> test_params(3);
  test_params[0] = params[0];
  test_params[1] = params[1]; 
  test_params[2] = params[2];
  
  tape.registerInput(test_params[0]);
  tape.newRecording();
  
  // Recreate system with registered params
  ActiveSys sys_test(test_params[0], test_params[1], test_params[2]);
  auto retrieved = sys_test.pars();
  
  ad_type result = retrieved[0] * 2.0;
  
  tape.registerOutput(result);
  xad::derivative(result) = 1.0;
  tape.computeAdjoints();
  
  double gradient = xad::derivative(test_params[0]);
  
  delete sys_active;
  
  return Rcpp::List::create(
    Rcpp::Named("passive_params") = params,
    Rcpp::Named("active_params") = params_active,
    Rcpp::Named("gradient_test") = gradient,
    Rcpp::Named("gradient_should_be") = 2.0
  );
}