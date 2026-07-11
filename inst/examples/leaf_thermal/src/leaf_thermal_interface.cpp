// [[Rcpp::depends(Rcpp, odelia)]]
// [[Rcpp::plugins(cpp20)]]

/* This file defines uses Rcpp to create an interface for the leaf thermal model, including
- parameter struct: LeafThermalPars
- model system: LeafThermalSystem
- drivers: drivers::Drivers
- ODE solver: ode::Solver<LeafThermalSystem>

Most are straightforward general wrappers around the underlying C++ methods and not specific to the leaf thermal model.

Throughout the functions below, we

- use functions names like Class_method to indicate method 'method' for class 'Class'. You should refer to these functions for details on underlying C++ methods.
- use Rcpp::XPtr to manage pointers to C++ objects created in R and passed back to C++.
- convert between Rcpp types (NumericVector, List, DataFrame) and C++ types (std::vector, structs) as needed.
- provide error checking for input sizes and validity where appropriate.
- return results in R-friendly formats.
- all functions are exported to R using the [[Rcpp::export]] attribute.

These functions are intended only as interface and are called via a corresponding R interface,to provide a user-friendly R interface. See `R/leaf_thermal_interface.R` for details.
*/

#include <odelia/solver_interface.hpp>
#include "leaf_thermal_system.hpp"

using namespace Rcpp;
using namespace odelia;

// Define types for Leaf Thermal system
typedef LeafThermalSystem<double> SystemType;
typedef LeafThermalSystem<xad::adj<double>::active_type> ActiveSystemType;

// Helper to get system pointer
inline Rcpp::XPtr<SystemType> get_LeafThermalSystem(SEXP xp) {
  return Rcpp::XPtr<SystemType>(xp);
}

// System interface (Leaf-specific)

inline List leaf_pars_to_list(const LeafThermalPars &p) {
  return List::create(
      _["k_H"] = p.k_H,
      _["g_tr_max"] = p.g_tr_max,
      _["m_tr"] = p.m_tr,
      _["T_tr_mid"] = p.T_tr_mid);
}

inline LeafThermalPars leaf_pars_from_list(const List &L) {
  LeafThermalPars p;

  if (L.containsElementNamed("k_H"))
    p.k_H = as<double>(L["k_H"]);
  if (L.containsElementNamed("g_tr_max"))
    p.g_tr_max = as<double>(L["g_tr_max"]);
  if (L.containsElementNamed("m_tr"))
    p.m_tr = as<double>(L["m_tr"]);
  if (L.containsElementNamed("T_tr_mid"))
    p.T_tr_mid = as<double>(L["T_tr_mid"]);

  return p;
}

// [[Rcpp::export]]
Rcpp::List LeafThermalSystemPars() {
  LeafThermalPars p;
  return leaf_pars_to_list(p);
}

// [[Rcpp::export]]
SEXP LeafThermalSystem_new(Rcpp::List pars_list, SEXP drivers_xp) {
  LeafThermalPars pars = leaf_pars_from_list(pars_list);
  auto drivers = get_Drivers(drivers_xp);
  Rcpp::XPtr<SystemType> ptr(new SystemType(pars, *drivers), true);
  return ptr;
}

// [[Rcpp::export]]
Rcpp::NumericVector LeafThermalSystem_pars(SEXP LeafThermalSystem_xp) {
  auto lor = get_LeafThermalSystem(LeafThermalSystem_xp);
  std::vector<double> p = lor->get_pars();
  return Rcpp::wrap(p);
}

// [[Rcpp::export]]
void LeafThermalSystem_initialize_drivers(SEXP LeafThermalSystem_xp, SEXP drivers_xp) {
  auto drv = get_Drivers(drivers_xp);
  auto lor = get_LeafThermalSystem(LeafThermalSystem_xp);
  lor->initialize_drivers(*drv);
}

// [[Rcpp::export]]
void LeafThermalSystem_set_state(SEXP LeafThermalSystem_xp, Rcpp::NumericVector y, double time) {
  auto lor = get_LeafThermalSystem(LeafThermalSystem_xp);
  if (y.size() != lor->ode_size()) {
    Rcpp::stop("State vector size mismatch");
  }
  std::vector<double> tmp(y.begin(), y.end());
  lor->set_ode_state(tmp.begin(), time);
}

// [[Rcpp::export]]
void LeafThermalSystem_set_initial_state(SEXP LeafThermalSystem_xp, Rcpp::NumericVector y, double t0 = 0.0) {
  auto lor = get_LeafThermalSystem(LeafThermalSystem_xp);
  if (y.size() != lor->ode_size()) {
    Rcpp::stop("State vector size mismatch");
  }
  std::vector<double> tmp(y.begin(), y.end());
  lor->set_initial_state(tmp.begin(), t0);
}

// [[Rcpp::export]]
void LeafThermalSystem_set_params(SEXP LeafThermalSystem_xp, Rcpp::NumericVector params) {
  auto lor = get_LeafThermalSystem(LeafThermalSystem_xp);
  std::vector<double> tmp(params.begin(), params.end());
  lor->set_params(tmp.begin());
  lor->compute_ode_rates();
}

// [[Rcpp::export]]
void LeafThermalSystem_reset(SEXP LeafThermalSystem_xp) {
  auto lor = get_LeafThermalSystem(LeafThermalSystem_xp);
  lor->reset();
}

// [[Rcpp::export]]
Rcpp::NumericVector LeafThermalSystem_state(SEXP LeafThermalSystem_xp) {
  auto lor = get_LeafThermalSystem(LeafThermalSystem_xp);
  std::vector<double> tmp(lor->ode_size());
  lor->ode_state(tmp.begin());
  return Rcpp::wrap(tmp);
}

// [[Rcpp::export]]
Rcpp::NumericVector LeafThermalSystem_rates(SEXP LeafThermalSystem_xp) {
  auto lor = get_LeafThermalSystem(LeafThermalSystem_xp);
  std::vector<double> tmp(lor->ode_size());
  lor->ode_rates(tmp.begin());
  return Rcpp::wrap(tmp);
}

// [[Rcpp::export]]
Rcpp::NumericVector LeafThermalSystem_get_current_drivers(SEXP LeafThermalSystem_xp) {
  auto lor = get_LeafThermalSystem(LeafThermalSystem_xp);
  return Rcpp::wrap(lor->get_current_drivers());
}

// Solver interface (Leaf-specific creation, generic operations)

// Build a double Solver for a LeafThermalSystem.
// [[Rcpp::export]]
SEXP LeafSolver_new(SEXP system_xp, SEXP control_xp, SEXP drivers_xp) {
  Rcpp::XPtr<SystemType> sys(system_xp);
  Rcpp::XPtr<ode::OdeControl> ctrl(control_xp);
  return Rcpp::XPtr<ode::Solver<SystemType>>(
      new ode::Solver<SystemType>(*sys, *ctrl), true);
}

// All other Solver functions call the generic (double-only) templates.

// [[Rcpp::export]]
void LeafSolver_reset(SEXP solver_xp) {
  odelia::solver::Solver_reset_impl<SystemType>(solver_xp);
}

// [[Rcpp::export]]
double LeafSolver_time(SEXP solver_xp) {
  return odelia::solver::Solver_time_impl<SystemType>(solver_xp);
}

// [[Rcpp::export]]
Rcpp::NumericVector LeafSolver_state(SEXP solver_xp) {
  return odelia::solver::Solver_state_impl<SystemType>(solver_xp);
}

// [[Rcpp::export]]
Rcpp::NumericVector LeafSolver_times(SEXP solver_xp) {
  return odelia::solver::Solver_times_impl<SystemType>(solver_xp);
}

// [[Rcpp::export]]
void LeafSolver_set_state(SEXP solver_xp, Rcpp::NumericVector y, double time) {
  odelia::solver::Solver_set_state_impl<SystemType>(solver_xp, y, time);
}

// [[Rcpp::export]]
void LeafSolver_advance_adaptive(SEXP solver_xp, Rcpp::NumericVector times) {
  odelia::solver::Solver_advance_adaptive_impl<SystemType>(solver_xp, times);
}

// [[Rcpp::export]]
void LeafSolver_advance_fixed(SEXP solver_xp, Rcpp::NumericVector times) {
  odelia::solver::Solver_advance_fixed_impl<SystemType>(solver_xp, times);
}

// [[Rcpp::export]]
void LeafSolver_advance_euler(SEXP solver_xp, Rcpp::NumericVector times) {
  odelia::solver::Solver_advance_euler_impl<SystemType>(solver_xp, times);
}

// [[Rcpp::export]]
void LeafSolver_step(SEXP solver_xp) {
  odelia::solver::Solver_step_impl<SystemType>(solver_xp);
}

// [[Rcpp::export]]
bool LeafSolver_get_collect(SEXP solver_xp) {
  return odelia::solver::Solver_get_collect_impl<SystemType>(solver_xp);
}

// [[Rcpp::export]]
void LeafSolver_set_collect(SEXP solver_xp, bool x) {
  odelia::solver::Solver_set_collect_impl<SystemType>(solver_xp, x);
}

// [[Rcpp::export]]
std::size_t LeafSolver_get_history_size(SEXP solver_xp) {
  return odelia::solver::Solver_get_history_size_impl<SystemType>(solver_xp);
}

// [[Rcpp::export]]
Rcpp::DataFrame LeafSolver_get_history_step(SEXP solver_xp, std::size_t i) {
  return odelia::solver::Solver_get_history_step_impl<SystemType>(solver_xp, i);
}

// [[Rcpp::export]]
Rcpp::List LeafSolver_get_history(SEXP solver_xp) {
  return odelia::solver::Solver_get_history_impl<SystemType>(solver_xp);
}

// Value + least-squares gradient on the double handle. Observations are passed per
// call and owned by the functional; the solver holds no calibration state.
// [[Rcpp::export]]
Rcpp::List LeafSolver_value_and_gradient(SEXP solver_xp,
                         Rcpp::NumericVector times,
                         Rcpp::NumericMatrix observations,
                         Rcpp::IntegerVector obs_indices,
                         Rcpp::Nullable<Rcpp::NumericVector> ic = R_NilValue,
                         Rcpp::Nullable<Rcpp::NumericVector> params = R_NilValue) {
  return odelia::solver::Solver_value_and_gradient_impl<SystemType, ActiveSystemType>(
    solver_xp, ic, params, times, observations, obs_indices
  );
}
