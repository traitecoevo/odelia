/* This file defines uses Rcpp to create an interface for the lorenz model, including

- model system: LorenzSystem
- ODE solver: ode::Solver<LorenzSystem>

Most are straightforward general wrappers around the underlying C++ methods and not specific to the lorenz model.

Throughout the functions below, we

- use functions names like Class_method to indicate method 'method' for class 'Class'. You should refer to these functions for details on underlying C++ methods.
- use Rcpp::XPtr to manage pointers to C++ objects created in R and passed back to C++.
- convert between Rcpp types (NumericVector, List, DataFrame) and C++ types (std::vector, structs) as needed.
- provide error checking for input sizes and validity where appropriate.
- return results in R-friendly formats.
- all functions are exported to R using the [[Rcpp::export]] attribute.

These functions are intended only as interface and are called via a corresponding R interface,to provide a user-friendly R interface. See `R/leaf_thermal_interface.R` for details.
*/
/* Rcpp interface for LorenzSystem with XAD automatic differentiation support */

#include <Rcpp.h>
#include <XAD/XAD.hpp>
#include <odelia/ode_solver.hpp>
#include <odelia/ode_fit.hpp>
#include <examples/lorenz_system.hpp>

using namespace Rcpp;
using namespace odelia;

// Define types for convenience
typedef LorenzSystem<double> SystemType;
typedef LorenzSystem<xad::adj<double>::active_type> ActiveSystemType;

// Helpers to convert pointers
template<typename T>
inline Rcpp::XPtr<ode::Solver<T>> get_solver(SEXP xp) {
  return Rcpp::XPtr<ode::Solver<T>>(xp);
}

inline Rcpp::XPtr<SystemType> get_system(SEXP xp) {
  return Rcpp::XPtr<SystemType>(xp);
}

//-------------------------------------------------------------------------
// Solver interface

// [[Rcpp::export]]
SEXP Solver_new(SEXP system_xp, SEXP control_xp, bool active = false) {
  Rcpp::XPtr<SystemType> sys(system_xp);
  Rcpp::XPtr<ode::OdeControl> ctrl(control_xp);
  
  if (active) {
    // Initialise new AD compatible system
    auto params = sys->pars();
    
    std::vector<double> initial_state(sys->ode_size());
    sys->ode_initial_state(initial_state.begin());
    auto t0 = sys->ode_t0();
    
    auto* sys_active = new ActiveSystemType(params[0], params[1], params[2]);
    sys_active->set_initial_state(initial_state.begin(), t0);
    
    auto* solver = new ode::Solver<ActiveSystemType>(*sys_active, *ctrl);
    
    return Rcpp::XPtr<ode::Solver<ActiveSystemType>>(solver, true);
  }else {
    // Create regular solver
    auto* sys_copy = new SystemType(*sys);
    auto* solver = new ode::Solver<SystemType>(*sys_copy, *ctrl);
    
    return Rcpp::XPtr<ode::Solver<SystemType>>(solver, true);
  }
}

// [[Rcpp::export]]
void Solver_reset(SEXP solver_xp, bool active = false) {
  if (active) {
    get_solver<ActiveSystemType>(solver_xp)->reset();
  } else {
    get_solver<SystemType>(solver_xp)->reset();
  }
}

// [[Rcpp::export]]
double Solver_time(SEXP solver_xp, bool active = false) {
  if (active) {
    return get_solver<ActiveSystemType>(solver_xp)->time();
  } else {
    return get_solver<SystemType>(solver_xp)->time();
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector Solver_state(SEXP solver_xp, bool active = false) {
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

// [[Rcpp::export]]
Rcpp::NumericVector Solver_times(SEXP solver_xp, bool active = false) {
  if (active) {
    return Rcpp::wrap(get_solver<ActiveSystemType>(solver_xp)->times());
  } else {
    return Rcpp::wrap(get_solver<SystemType>(solver_xp)->times());
  }
}

// [[Rcpp::export]]
void Solver_set_state(SEXP solver_xp, Rcpp::NumericVector y, double time, bool active = false) {
  std::vector<double> yy(y.begin(), y.end());
  if (active) {
    get_solver<ActiveSystemType>(solver_xp)->set_state(yy, time);
  } else {
    get_solver<SystemType>(solver_xp)->set_state(yy, time);
  }
}

// [[Rcpp::export]]
void Solver_advance_adaptive(SEXP solver_xp, Rcpp::NumericVector times, bool active = false) {
  if (active) {
    Rcpp::stop("advance_adaptive() not supported for AD solvers. Use advance_fixed() with pre-computed schedule.");
  }
  std::vector<double> ts(times.begin(), times.end());
  get_solver<SystemType>(solver_xp)->advance_adaptive(ts);
}

// [[Rcpp::export]]
void Solver_advance_fixed(SEXP solver_xp, Rcpp::NumericVector times, bool active = false) {
  std::vector<double> ts(times.begin(), times.end());
  if (active) {
    get_solver<ActiveSystemType>(solver_xp)->advance_fixed(ts);
  } else {
    get_solver<SystemType>(solver_xp)->advance_fixed(ts);
  }
}

// [[Rcpp::export]]
void Solver_step(SEXP solver_xp, bool active = false) {
  if (active) {
    Rcpp::stop("step() not supported for AD solvers. Use advance_fixed() instead.");
  }
  get_solver<SystemType>(solver_xp)->step();
}

// [[Rcpp::export]]
bool Solver_get_collect(SEXP solver_xp, bool active = false) {
  if (active) {
    return get_solver<ActiveSystemType>(solver_xp)->get_collect();
  } else {
    return get_solver<SystemType>(solver_xp)->get_collect();
  }
}

// [[Rcpp::export]]
void Solver_set_collect(SEXP solver_xp, bool x, bool active = false) {
  if (active) {
    get_solver<ActiveSystemType>(solver_xp)->set_collect(x);
  } else {
    get_solver<SystemType>(solver_xp)->set_collect(x);
  }
}

// [[Rcpp::export]]
std::size_t Solver_get_history_size(SEXP solver_xp, bool active = false) {
  if (active) {
    return get_solver<ActiveSystemType>(solver_xp)->get_history_size();
  } else {
    return get_solver<SystemType>(solver_xp)->get_history_size();
  }
}

// Helper to get column names
CharacterVector get_column_names() {
  return CharacterVector::create("time", "x", "y", "z", "dxdt", "dydt", "dzdt");
}

// [[Rcpp::export]]
Rcpp::DataFrame Solver_get_history_step(SEXP solver_xp, std::size_t i, bool active = false) {
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
  
  CharacterVector names = get_column_names();
  Rcpp::List df_list(names.size());
  for (size_t j = 0; j < names.size(); ++j) {
    df_list[j] = out[j];
  }
  df_list.attr("names") = names;
  
  return Rcpp::DataFrame(df_list);
}

// [[Rcpp::export]]
Rcpp::List Solver_get_history(SEXP solver_xp, bool active = false) {
  int nrows;
  std::vector<std::vector<double>> cols;
  CharacterVector names = get_column_names();
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

//-------------------------------------------------------------------------
// System interface

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
// OdeControl interface

// [[Rcpp::export]]
SEXP OdeControl_new() {
  return Rcpp::XPtr<ode::OdeControl>(new ode::OdeControl(), true);
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
// Fitting interface

// [[Rcpp::export]]
void Solver_set_target(SEXP solver_xp, 
                      Rcpp::NumericVector times,
                      Rcpp::NumericMatrix target,
                      Rcpp::IntegerVector obs_indices,
                      bool active = false) {
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

// COMBINED FIT FUNCTION - replaces Solver_fit_ic and Solver_fit_params
// [[Rcpp::export]]
Rcpp::List Solver_fit(SEXP solver_xp,
                     Rcpp::Nullable<Rcpp::NumericVector> ic = R_NilValue,
                     Rcpp::Nullable<Rcpp::NumericVector> params = R_NilValue) {
    
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
    
    // Compute gradient using new unified function
    auto [loss, gradient] = ode::compute_gradient(*solver, ic_opt, params_opt);
    
    return Rcpp::List::create(
        Rcpp::Named("loss") = loss,
        Rcpp::Named("gradient") = Rcpp::wrap(gradient)
    );
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
