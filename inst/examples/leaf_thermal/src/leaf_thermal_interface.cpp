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

#include <Rcpp.h>
#include <XAD/XAD.hpp>
#include <odelia/ode_solver.hpp>
#include <odelia/ode_fit.hpp>
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

// Solver creation - Leaf-specific (must know LeafThermalSystem type)
// [[Rcpp::export]]
SEXP LeafSolver_new(SEXP system_xp, SEXP control_xp, SEXP drivers_xp, bool active = false) {
  Rcpp::XPtr<SystemType> sys(system_xp);
  Rcpp::XPtr<ode::OdeControl> ctrl(control_xp);
  Rcpp::XPtr<drivers::Drivers> drv(drivers_xp);

  if (active) {
    auto pars = sys->get_pars();
    
    std::vector<double> initial_state(sys->ode_size());
    sys->ode_initial_state(initial_state.begin());
    auto t0 = sys->ode_t0();
    
    auto* sys_active = new ActiveSystemType(LeafThermalPars{pars[0], pars[1], pars[2], pars[3]}, *drv);
    sys_active->set_initial_state(initial_state.begin(), t0);
    
    auto* solver = new ode::Solver<ActiveSystemType>(*sys_active, *ctrl);
    return Rcpp::XPtr<ode::Solver<ActiveSystemType>>(solver, true);
  } else {
    auto* sys_copy = new SystemType(*sys);
    auto* solver = new ode::Solver<SystemType>(*sys_copy, *ctrl);
    return Rcpp::XPtr<ode::Solver<SystemType>>(solver, true);
  }
}

// Helper to get column names (Leaf-specific)
inline Rcpp::CharacterVector get_leaf_column_names(SEXP solver_xp, bool active) {
  if (active) {
    auto solver = odelia::solver::get_solver<ActiveSystemType>(solver_xp);
    return Rcpp::wrap(solver->get_system().record_colnames());
  } else {
    auto solver = odelia::solver::get_solver<SystemType>(solver_xp);
    return Rcpp::wrap(solver->get_system().record_colnames());
  }
}

// All other Solver functions call generic templates with LeafThermalSystem types

// [[Rcpp::export]]
void LeafSolver_reset(SEXP solver_xp, bool active = false) {
  odelia::solver::Solver_reset_impl<SystemType, ActiveSystemType>(solver_xp, active);
}

// [[Rcpp::export]]
double LeafSolver_time(SEXP solver_xp, bool active = false) {
  return odelia::solver::Solver_time_impl<SystemType, ActiveSystemType>(solver_xp, active);
}

// [[Rcpp::export]]
Rcpp::NumericVector LeafSolver_state(SEXP solver_xp, bool active = false) {
  return odelia::solver::Solver_state_impl<SystemType, ActiveSystemType>(solver_xp, active);
}

// [[Rcpp::export]]
Rcpp::NumericVector LeafSolver_times(SEXP solver_xp, bool active = false) {
  return odelia::solver::Solver_times_impl<SystemType, ActiveSystemType>(solver_xp, active);
}

// [[Rcpp::export]]
void LeafSolver_set_state(SEXP solver_xp, Rcpp::NumericVector y, double time, bool active = false) {
  odelia::solver::Solver_set_state_impl<SystemType, ActiveSystemType>(solver_xp, y, time, active);
}

// [[Rcpp::export]]
void LeafSolver_advance_adaptive(SEXP solver_xp, Rcpp::NumericVector times, bool active = false) {
  odelia::solver::Solver_advance_adaptive_impl<SystemType, ActiveSystemType>(solver_xp, times, active);
}

// [[Rcpp::export]]
void LeafSolver_advance_fixed(SEXP solver_xp, Rcpp::NumericVector times, bool active = false) {
  odelia::solver::Solver_advance_fixed_impl<SystemType, ActiveSystemType>(solver_xp, times, active);
}

// [[Rcpp::export]]
void LeafSolver_step(SEXP solver_xp, bool active = false) {
  odelia::solver::Solver_step_impl<SystemType, ActiveSystemType>(solver_xp, active);
}

// [[Rcpp::export]]
bool LeafSolver_get_collect(SEXP solver_xp, bool active = false) {
  return odelia::solver::Solver_get_collect_impl<SystemType, ActiveSystemType>(solver_xp, active);
}

// [[Rcpp::export]]
void LeafSolver_set_collect(SEXP solver_xp, bool x, bool active = false) {
  odelia::solver::Solver_set_collect_impl<SystemType, ActiveSystemType>(solver_xp, x, active);
}

// [[Rcpp::export]]
std::size_t LeafSolver_get_history_size(SEXP solver_xp, bool active = false) {
  return odelia::solver::Solver_get_history_size_impl<SystemType, ActiveSystemType>(solver_xp, active);
}

// [[Rcpp::export]]
Rcpp::DataFrame LeafSolver_get_history_step(SEXP solver_xp, std::size_t i, bool active = false) {
  return odelia::solver::Solver_get_history_step_impl<SystemType, ActiveSystemType>(
    solver_xp, i, get_leaf_column_names(solver_xp, active), active
  );
}

// [[Rcpp::export]]
Rcpp::List LeafSolver_get_history(SEXP solver_xp, bool active = false) {
  return odelia::solver::Solver_get_history_impl<SystemType, ActiveSystemType>(
    solver_xp, get_leaf_column_names(solver_xp, active), active
  );
}

// [[Rcpp::export]]
void LeafSolver_set_target(SEXP solver_xp, 
                          Rcpp::NumericVector times,
                          Rcpp::NumericMatrix target,
                          Rcpp::IntegerVector obs_indices,
                          bool active = false) {
  odelia::solver::Solver_set_target_impl<SystemType, ActiveSystemType>(
    solver_xp, times, target, obs_indices, active
  );
}

// [[Rcpp::export]]
Rcpp::List LeafSolver_fit(SEXP solver_xp,
                         Rcpp::Nullable<Rcpp::NumericVector> ic = R_NilValue,
                         Rcpp::Nullable<Rcpp::NumericVector> params = R_NilValue) {
  return odelia::solver::Solver_fit_impl<SystemType, ActiveSystemType>(
    solver_xp, ic, params
  );
}
