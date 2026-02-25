// [[Rcpp::depends(Rcpp, odelia)]]
// [[Rcpp::plugins(cpp20)]]

#include <Rcpp.h>
#include <XAD/XAD.hpp>
#include <odelia/ode_solver.hpp>
#include <odelia/ode_fit.hpp>
#include "leaf_thermal_system.hpp"

using namespace Rcpp;
using namespace odelia;

// Define types for convenience
typedef LeafThermalSystem<double> SystemType;
typedef LeafThermalSystem<xad::adj<double>::active_type> ActiveSystemType;

// Helpers to convert pointers from SEXP to XPtr
template<typename T>
inline Rcpp::XPtr<ode::Solver<T>> get_leaf_solver(SEXP xp) {
  return Rcpp::XPtr<ode::Solver<T>>(xp);
}

inline Rcpp::XPtr<SystemType> get_LeafThermalSystem(SEXP xp) {
  return Rcpp::XPtr<SystemType>(xp);
}

inline Rcpp::XPtr<drivers::Drivers> get_Drivers(SEXP xp) {
  return Rcpp::XPtr<drivers::Drivers>(xp);
}

inline Rcpp::XPtr<ode::OdeControl> get_OdeControl(SEXP xp){
  return Rcpp::XPtr<ode::OdeControl>(xp);
}

// Rcpp interface for ode::Solver<LeafThermalSystem>

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

// [[Rcpp::export]]
void LeafSolver_reset(SEXP solver_xp, bool active = false) {
  if (active) {
    get_leaf_solver<ActiveSystemType>(solver_xp)->reset();
  } else {
    get_leaf_solver<SystemType>(solver_xp)->reset();
  }
}

// [[Rcpp::export]]
double LeafSolver_time(SEXP solver_xp, bool active = false) {
  if (active) {
    return get_leaf_solver<ActiveSystemType>(solver_xp)->time();
  } else {
    return get_leaf_solver<SystemType>(solver_xp)->time();
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector LeafSolver_state(SEXP solver_xp, bool active = false) {
  if (active) {
    auto state = get_leaf_solver<ActiveSystemType>(solver_xp)->state();
    std::vector<double> result(state.size());
    for (size_t i = 0; i < state.size(); ++i) {
      result[i] = xad::value(state[i]);
    }
    return Rcpp::wrap(result);
  } else {
    auto state = get_leaf_solver<SystemType>(solver_xp)->state();
    return Rcpp::wrap(state);
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector LeafSolver_times(SEXP solver_xp, bool active = false) {
  if (active) {
    return Rcpp::wrap(get_leaf_solver<ActiveSystemType>(solver_xp)->times());
  } else {
    return Rcpp::wrap(get_leaf_solver<SystemType>(solver_xp)->times());
  }
}

// [[Rcpp::export]]
void LeafSolver_set_state(SEXP solver_xp, Rcpp::NumericVector y, double time, bool active = false) {
  std::vector<double> yy(y.begin(), y.end());
  if (active) {
    get_leaf_solver<ActiveSystemType>(solver_xp)->set_state(yy, time);
  } else {
    get_leaf_solver<SystemType>(solver_xp)->set_state(yy, time);
  }
}

// [[Rcpp::export]]
void LeafSolver_advance_adaptive(SEXP solver_xp, Rcpp::NumericVector times, bool active = false) {
  if (active) {
    Rcpp::stop("advance_adaptive() not supported for AD solvers. Use advance_fixed() with pre-computed schedule.");
  }
  std::vector<double> ts(times.begin(), times.end());
  get_leaf_solver<SystemType>(solver_xp)->advance_adaptive(ts);
}

// [[Rcpp::export]]
void LeafSolver_advance_fixed(SEXP solver_xp, Rcpp::NumericVector times, bool active = false) {
  std::vector<double> ts(times.begin(), times.end());
  if (active) {
    get_leaf_solver<ActiveSystemType>(solver_xp)->advance_fixed(ts);
  } else {
    get_leaf_solver<SystemType>(solver_xp)->advance_fixed(ts);
  }
}

// [[Rcpp::export]]
void LeafSolver_step(SEXP solver_xp, bool active = false) {
  if (active) {
    Rcpp::stop("step() not supported for AD solvers. Use advance_fixed() instead.");
  }
  get_leaf_solver<SystemType>(solver_xp)->step();
}

// [[Rcpp::export]]
bool LeafSolver_get_collect(SEXP solver_xp, bool active = false) {
  if (active) {
    return get_leaf_solver<ActiveSystemType>(solver_xp)->get_collect();
  } else {
    return get_leaf_solver<SystemType>(solver_xp)->get_collect();
  }
}

// [[Rcpp::export]]
void LeafSolver_set_collect(SEXP solver_xp, bool x, bool active = false) {
  if (active) {
    get_leaf_solver<ActiveSystemType>(solver_xp)->set_collect(x);
  } else {
    get_leaf_solver<SystemType>(solver_xp)->set_collect(x);
  }
}

// [[Rcpp::export]]
std::size_t LeafSolver_get_history_size(SEXP solver_xp, bool active = false) {
  if (active) {
    return get_leaf_solver<ActiveSystemType>(solver_xp)->get_history_size();
  } else {
    return get_leaf_solver<SystemType>(solver_xp)->get_history_size();
  }
}

// [[Rcpp::export]]
Rcpp::DataFrame LeafSolver_get_history_step(SEXP solver_xp, std::size_t i, bool active = false) {
  std::vector<double> out;
  std::vector<std::string> names;
  
  if (active) {
    auto solver = get_leaf_solver<ActiveSystemType>(solver_xp);
    if (i >= solver->get_history_size()) {
      Rcpp::stop("Index out of bounds");
    }
    out = solver->get_history_step(i).record_step();
    names = solver->get_system().record_colnames();
  } else {
    auto solver = get_leaf_solver<SystemType>(solver_xp);
    if (i >= solver->get_history_size()) {
      Rcpp::stop("Index out of bounds");
    }
    out = solver->get_history_step(i).record_step();
    names = solver->get_system().record_colnames();
  }

  Rcpp::List df_list(names.size());
  for (size_t j = 0; j < names.size(); ++j) {
    df_list[j] = Rcpp::wrap(out[j]);
  }
  df_list.attr("names") = names;

  return Rcpp::DataFrame(df_list);
}

// [[Rcpp::export]]
Rcpp::List LeafSolver_get_history(SEXP solver_xp, bool active = false) {
  int nrows;
  std::vector<std::vector<double>> cols;
  std::vector<std::string> names;
  int ncols;

  if (active) {
    auto solver = get_leaf_solver<ActiveSystemType>(solver_xp);
    nrows = solver->get_history_size();
    names = solver->get_system().record_colnames();
    ncols = names.size();
    cols.resize(ncols);
    for (auto& col : cols) col.reserve(nrows);
    
    for (size_t i = 0; i < nrows; ++i) {
      auto row = solver->get_history_step(i).record_step();
      for (size_t j = 0; j < ncols; ++j) {
        cols[j].push_back(row[j]);
      }
    }
  } else {
    auto solver = get_leaf_solver<SystemType>(solver_xp);
    nrows = solver->get_history_size();
    names = solver->get_system().record_colnames();
    ncols = names.size();
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

// Fitting interface

// [[Rcpp::export]]
void LeafSolver_set_target(SEXP solver_xp, 
                      Rcpp::NumericVector times,
                      Rcpp::NumericMatrix target,
                      Rcpp::IntegerVector obs_indices,
                      bool active = false) {
    std::vector<double> times_vec(times.begin(), times.end());
    
    int nrows = target.nrow(), ncols = target.ncol();
    std::vector<std::vector<double>> targets_vec(nrows);
    for (int i = 0; i < nrows; ++i) {
        targets_vec[i].resize(ncols);
        for (int j = 0; j < ncols; ++j) {
            targets_vec[i][j] = target(i, j);
        }
    }
    
    std::vector<size_t> obs_idx_vec(obs_indices.size());
    for (size_t i = 0; i < obs_indices.size(); ++i) {
        obs_idx_vec[i] = obs_indices[i] - 1;
    }
    
    if (active) {
        get_leaf_solver<ActiveSystemType>(solver_xp)->set_target(times_vec, targets_vec, obs_idx_vec);
    } else {
        get_leaf_solver<SystemType>(solver_xp)->set_target(times_vec, targets_vec, obs_idx_vec);
    }
}

// [[Rcpp::export]]
Rcpp::List LeafSolver_fit(SEXP solver_xp,
                     Rcpp::Nullable<Rcpp::NumericVector> ic = R_NilValue,
                     Rcpp::Nullable<Rcpp::NumericVector> params = R_NilValue) {
    
    if (ic.isNull() && params.isNull()) {
        Rcpp::stop("Must provide at least one of 'ic' or 'params'");
    }
    
    auto solver = get_leaf_solver<ActiveSystemType>(solver_xp);
    
    if (solver->fit_times().empty()) {
        Rcpp::stop("Must call set_target() before LeafSolver_fit()");
    }
    
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
    
    auto [loss, gradient] = ode::compute_gradient(*solver, ic_opt, params_opt);
    
    return Rcpp::List::create(
        Rcpp::Named("loss") = loss,
        Rcpp::Named("gradient") = Rcpp::wrap(gradient)
    );
}

// Rcpp interface for LeafThermalSystem

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
