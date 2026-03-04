/* Generic ODE interface for odelia package
 * 
 * This file provides generic functions that can be used with any ODE system.
 * It includes:
 * - Drivers: External forcing data (constant or time-varying)
 * - OdeControl: Solver settings and tolerances
 * - Generic Solver functions: Templated implementations that work with any System type
 * 
 * These functions are system-agnostic and can be used by any model
 * (Lorenz, Leaf Thermal, ATLAS, etc.)
 * 
 * System-specific interfaces (lorenz_interface.cpp, leaf_thermal_interface.cpp)
 * call these templated functions with their specific System types.
 */

#include <Rcpp.h>
#include <XAD/XAD.hpp>
#include <odelia/drivers.hpp>
#include <odelia/ode_control.hpp>
#include <odelia/ode_solver.hpp>
#include <odelia/ode_fit.hpp>

using namespace Rcpp;
using namespace odelia;

//-------------------------------------------------------------------------
// OdeControl interface

// Helper to get OdeControl pointer
inline Rcpp::XPtr<ode::OdeControl> get_OdeControl(SEXP xp) {
  return Rcpp::XPtr<ode::OdeControl>(xp);
}

// [[Rcpp::export]]
SEXP OdeControl_new() {
  return Rcpp::XPtr<ode::OdeControl>(new ode::OdeControl(), true);
}

// [[Rcpp::export]]
Rcpp::List OdeControl_get_controls(SEXP control_xp) {
  auto ctrl = get_OdeControl(control_xp);
  return Rcpp::List::create(
    Rcpp::Named("tol_abs") = ctrl->tol_abs,
    Rcpp::Named("tol_rel") = ctrl->tol_rel,
    Rcpp::Named("a_y") = ctrl->a_y,
    Rcpp::Named("a_dydt") = ctrl->a_dydt,
    Rcpp::Named("step_size_min") = ctrl->step_size_min,
    Rcpp::Named("step_size_max") = ctrl->step_size_max,
    Rcpp::Named("step_size_initial") = ctrl->step_size_initial
  );
}

// [[Rcpp::export]]
void OdeControl_set_controls(SEXP control_xp,
                             double tol_abs, double tol_rel,
                             double a_y, double a_dydt,
                             double step_size_min, double step_size_max,
                             double step_size_initial) {
  auto ctrl = get_OdeControl(control_xp);
  ctrl->tol_abs = tol_abs;
  ctrl->tol_rel = tol_rel;
  ctrl->a_y = a_y;
  ctrl->a_dydt = a_dydt;
  ctrl->step_size_min = step_size_min;
  ctrl->step_size_max = step_size_max;
  ctrl->step_size_initial = step_size_initial;
}

// [[Rcpp::export]]
void OdeControl_set_tol_abs(SEXP control_xp, double tol_abs) {
  get_OdeControl(control_xp)->tol_abs = tol_abs;
}

// [[Rcpp::export]]
void OdeControl_set_tol_rel(SEXP control_xp, double tol_rel) {
  get_OdeControl(control_xp)->tol_rel = tol_rel;
}

// [[Rcpp::export]]
void OdeControl_set_a_y(SEXP control_xp, double a_y) {
  get_OdeControl(control_xp)->a_y = a_y;
}

// [[Rcpp::export]]
void OdeControl_set_a_dydt(SEXP control_xp, double a_dydt) {
  get_OdeControl(control_xp)->a_dydt = a_dydt;
}

// [[Rcpp::export]]
void OdeControl_set_step_size_min(SEXP control_xp, double step_size_min) {
  get_OdeControl(control_xp)->step_size_min = step_size_min;
}

// [[Rcpp::export]]
void OdeControl_set_step_size_max(SEXP control_xp, double step_size_max) {
  get_OdeControl(control_xp)->step_size_max = step_size_max;
}

// [[Rcpp::export]]
void OdeControl_set_step_size_initial(SEXP control_xp, double step_size_initial) {
  get_OdeControl(control_xp)->step_size_initial = step_size_initial;
}

//-------------------------------------------------------------------------
// Drivers interface

// Helper to get Drivers pointer
inline Rcpp::XPtr<drivers::Drivers> get_Drivers(SEXP xp) {
  return Rcpp::XPtr<drivers::Drivers>(xp);
}

// [[Rcpp::export]]
SEXP Drivers_new() {
  Rcpp::XPtr<drivers::Drivers> ptr(new drivers::Drivers(), true);
  return ptr;
}

// [[Rcpp::export]]
void Drivers_set_constant(SEXP drivers_xp, std::string driver_name, double k) {
  auto drv = get_Drivers(drivers_xp);
  drv->set_constant(driver_name, k);
}

// [[Rcpp::export]]
void Drivers_set_variable(SEXP drivers_xp, std::string driver_name,
                          Rcpp::NumericVector x, Rcpp::NumericVector y) {
  auto drv = get_Drivers(drivers_xp);
  std::vector<double> x_vec(x.begin(), x.end());
  std::vector<double> y_vec(y.begin(), y.end());
  drv->set_variable(driver_name, x_vec, y_vec);
}

// [[Rcpp::export]]
void Drivers_set_extrapolate(SEXP drivers_xp, std::string driver_name, bool extrapolate) {
  auto drv = get_Drivers(drivers_xp);
  drv->set_extrapolate(driver_name, extrapolate);
}

// [[Rcpp::export]]
double Drivers_evaluate(SEXP drivers_xp, std::string driver_name, double x) {
  auto drv = get_Drivers(drivers_xp);
  return drv->evaluate(driver_name, x);
}

// [[Rcpp::export]]
Rcpp::NumericVector Drivers_evaluate_range(SEXP drivers_xp, std::string driver_name,
                                           Rcpp::NumericVector x) {
  auto drv = get_Drivers(drivers_xp);
  std::vector<double> x_vec(x.begin(), x.end());
  std::vector<double> result = drv->evaluate_range(driver_name, x_vec);
  return Rcpp::wrap(result);
}

// [[Rcpp::export]]
Rcpp::CharacterVector Drivers_get_names(SEXP drivers_xp) {
  auto drv = get_Drivers(drivers_xp);
  std::vector<std::string> names = drv->get_names();
  return Rcpp::wrap(names);
}

// [[Rcpp::export]]
void Drivers_clear(SEXP drivers_xp) {
  auto drv = get_Drivers(drivers_xp);
  drv->clear();
}

//-------------------------------------------------------------------------
// Generic Solver interface (templated implementations)
//
// These functions work with any System type. System-specific interfaces
// call these with their specific types (e.g., LorenzSystem, LeafThermalSystem)

namespace odelia {
namespace solver {

// Helper to get solver pointer (templated)
template<typename T>
inline Rcpp::XPtr<ode::Solver<T>> get_solver(SEXP xp) {
  return Rcpp::XPtr<ode::Solver<T>>(xp);
}

// Generic Solver_reset
template<typename SystemType, typename ActiveSystemType>
void Solver_reset_impl(SEXP solver_xp, bool active) {
  if (active) {
    get_solver<ActiveSystemType>(solver_xp)->reset();
  } else {
    get_solver<SystemType>(solver_xp)->reset();
  }
}

// Generic Solver_time
template<typename SystemType, typename ActiveSystemType>
double Solver_time_impl(SEXP solver_xp, bool active) {
  if (active) {
    return get_solver<ActiveSystemType>(solver_xp)->time();
  } else {
    return get_solver<SystemType>(solver_xp)->time();
  }
}

// Generic Solver_state
template<typename SystemType, typename ActiveSystemType>
Rcpp::NumericVector Solver_state_impl(SEXP solver_xp, bool active) {
  if (active) {
    auto state = get_solver<ActiveSystemType>(solver_xp)->state();
    std::vector<double> result(state.size());
    for (size_t i = 0; i < state.size(); ++i) {
      result[i] = xad::value(state[i]);
    }
    return Rcpp::wrap(result);
  } else {
    auto state = get_solver<SystemType>(solver_xp)->state();
    return Rcpp::wrap(state);
  }
}

// Generic Solver_times
template<typename SystemType, typename ActiveSystemType>
Rcpp::NumericVector Solver_times_impl(SEXP solver_xp, bool active) {
  if (active) {
    return Rcpp::wrap(get_solver<ActiveSystemType>(solver_xp)->times());
  } else {
    return Rcpp::wrap(get_solver<SystemType>(solver_xp)->times());
  }
}

// Generic Solver_set_state
template<typename SystemType, typename ActiveSystemType>
void Solver_set_state_impl(SEXP solver_xp, Rcpp::NumericVector y, double time, bool active) {
  std::vector<double> yy(y.begin(), y.end());
  if (active) {
    get_solver<ActiveSystemType>(solver_xp)->set_state(yy, time);
  } else {
    get_solver<SystemType>(solver_xp)->set_state(yy, time);
  }
}

// Generic Solver_advance_adaptive
template<typename SystemType, typename ActiveSystemType>
void Solver_advance_adaptive_impl(SEXP solver_xp, Rcpp::NumericVector times, bool active) {
  if (active) {
    Rcpp::stop("advance_adaptive() not supported for AD solvers. Use advance_fixed() with pre-computed schedule.");
  }
  std::vector<double> ts(times.begin(), times.end());
  get_solver<SystemType>(solver_xp)->advance_adaptive(ts);
}

// Generic Solver_advance_fixed
template<typename SystemType, typename ActiveSystemType>
void Solver_advance_fixed_impl(SEXP solver_xp, Rcpp::NumericVector times, bool active) {
  std::vector<double> ts(times.begin(), times.end());
  if (active) {
    get_solver<ActiveSystemType>(solver_xp)->advance_fixed(ts);
  } else {
    get_solver<SystemType>(solver_xp)->advance_fixed(ts);
  }
}

// Generic Solver_step
template<typename SystemType, typename ActiveSystemType>
void Solver_step_impl(SEXP solver_xp, bool active) {
  if (active) {
    Rcpp::stop("step() not supported for AD solvers. Use advance_fixed() instead.");
  }
  get_solver<SystemType>(solver_xp)->step();
}

// Generic Solver_get_collect
template<typename SystemType, typename ActiveSystemType>
bool Solver_get_collect_impl(SEXP solver_xp, bool active) {
  if (active) {
    return get_solver<ActiveSystemType>(solver_xp)->get_collect();
  } else {
    return get_solver<SystemType>(solver_xp)->get_collect();
  }
}

// Generic Solver_set_collect
template<typename SystemType, typename ActiveSystemType>
void Solver_set_collect_impl(SEXP solver_xp, bool x, bool active) {
  if (active) {
    get_solver<ActiveSystemType>(solver_xp)->set_collect(x);
  } else {
    get_solver<SystemType>(solver_xp)->set_collect(x);
  }
}

// Generic Solver_get_history_size
template<typename SystemType, typename ActiveSystemType>
std::size_t Solver_get_history_size_impl(SEXP solver_xp, bool active) {
  if (active) {
    return get_solver<ActiveSystemType>(solver_xp)->get_history_size();
  } else {
    return get_solver<SystemType>(solver_xp)->get_history_size();
  }
}

// Generic Solver_get_history_step
template<typename SystemType, typename ActiveSystemType>
Rcpp::DataFrame Solver_get_history_step_impl(SEXP solver_xp, std::size_t i, 
                                              CharacterVector names, bool active) {
  std::vector<double> out;
  
  if (active) {
    auto solver = get_solver<ActiveSystemType>(solver_xp);
    if (i >= solver->get_history_size()) {
      Rcpp::stop("Index out of bounds");
    }
    out = solver->get_history_step(i).record_step();
  } else {
    auto solver = get_solver<SystemType>(solver_xp);
    if (i >= solver->get_history_size()) {
      Rcpp::stop("Index out of bounds");
    }
    out = solver->get_history_step(i).record_step();
  }
  
  Rcpp::List df_list(names.size());
  for (size_t j = 0; j < names.size(); ++j) {
    df_list[j] = out[j];
  }
  df_list.attr("names") = names;
  
  return Rcpp::DataFrame(df_list);
}

// Generic Solver_get_history
template<typename SystemType, typename ActiveSystemType>
Rcpp::List Solver_get_history_impl(SEXP solver_xp, CharacterVector names, bool active) {
  int nrows;
  std::vector<std::vector<double>> cols;
  int ncols = names.size();
  
  if (active) {
    auto solver = get_solver<ActiveSystemType>(solver_xp);
    nrows = solver->get_history_size();
    cols.resize(ncols);
    for (auto& col : cols) col.reserve(nrows);
    
    for (size_t i = 0; i < nrows; ++i) {
      auto row = solver->get_history_step(i).record_step();
      for (size_t j = 0; j < ncols; ++j) {
        cols[j].push_back(row[j]);
      }
    }
  } else {
    auto solver = get_solver<SystemType>(solver_xp);
    nrows = solver->get_history_size();
    cols.resize(ncols);
    for (auto& col : cols) col.reserve(nrows);
    
    for (size_t i = 0; i < nrows; ++i) {
      auto row = solver->get_history_step(i).record_step();
      for (size_t j = 0; j < ncols; ++j) {
        cols[j].push_back(row[j]);
      }
    }
  }
  
  Rcpp::List out(ncols);
  for (size_t j = 0; j < ncols; ++j) {
    out[j] = NumericVector(cols[j].begin(), cols[j].end());
  }
  out.attr("names") = names;
  
  return DataFrame(out);
}

// Generic Solver_set_target
template<typename SystemType, typename ActiveSystemType>
void Solver_set_target_impl(SEXP solver_xp, 
                            Rcpp::NumericVector times,
                            Rcpp::NumericMatrix target,
                            Rcpp::IntegerVector obs_indices,
                            bool active) {
  // Convert times
  std::vector<double> times_vec(times.begin(), times.end());
  
  // Convert matrix
  int nrows = target.nrow(), ncols = target.ncol();
  std::vector<std::vector<double>> targets_vec(nrows);
  for (int i = 0; i < nrows; ++i) {
    targets_vec[i].resize(ncols);
    for (int j = 0; j < ncols; ++j) {
      targets_vec[i][j] = target(i, j);
    }
  }
  
  // Convert to 0-based indices
  std::vector<size_t> obs_idx_vec(obs_indices.size());
  for (size_t i = 0; i < obs_indices.size(); ++i) {
    obs_idx_vec[i] = obs_indices[i] - 1;
  }
  
  if (active) {
    get_solver<ActiveSystemType>(solver_xp)->set_target(times_vec, targets_vec, obs_idx_vec);
  } else {
    get_solver<SystemType>(solver_xp)->set_target(times_vec, targets_vec, obs_idx_vec);
  }
}

// Generic Solver_fit
template<typename SystemType, typename ActiveSystemType>
Rcpp::List Solver_fit_impl(SEXP solver_xp,
                           Rcpp::Nullable<Rcpp::NumericVector> ic,
                           Rcpp::Nullable<Rcpp::NumericVector> params) {
  // At least one must be provided
  if (ic.isNull() && params.isNull()) {
    Rcpp::stop("Must provide at least one of 'ic' or 'params'");
  }
  
  // Get AD solver (fit always uses AD)
  auto solver = get_solver<ActiveSystemType>(solver_xp);
  
  if (solver->fit_times().empty()) {
    Rcpp::stop("Must call set_target() before Solver_fit()");
  }
  
  // Convert to std::optional
  std::optional<std::vector<double>> ic_opt;
  if (!ic.isNull()) {
    Rcpp::NumericVector ic_vec(ic);
    ic_opt = std::vector<double>(ic_vec.begin(), ic_vec.end());
  }
  
  std::optional<std::vector<double>> params_opt;
  if (!params.isNull()) {
    Rcpp::NumericVector params_vec(params);
    params_opt = std::vector<double>(params_vec.begin(), params_vec.end());
  }
  
  // Compute gradient using unified function
  auto [loss, gradient] = ode::compute_gradient(*solver, ic_opt, params_opt);
  
  return Rcpp::List::create(
    Rcpp::Named("loss") = loss,
    Rcpp::Named("gradient") = Rcpp::wrap(gradient)
  );
}

} // namespace solver
} // namespace odelia