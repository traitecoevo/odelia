/* Generic Solver interface templates for odelia package
 * 
 * This header provides templated implementations of generic Solver functions
 * that work with any System type. System-specific interfaces include this
 * header and instantiate the templates with their specific types.
 * 
 * These templates are defined inline in the header to avoid linking issues.
 */

#ifndef ODELIA_SOLVER_INTERFACE_HPP_
#define ODELIA_SOLVER_INTERFACE_HPP_

#include <Rcpp.h>
#include <XAD/XAD.hpp>
#include <odelia/ode_solver.hpp>
#include <odelia/ode_fit.hpp>

namespace odelia {
namespace solver {

// Helper to get solver pointer (templated)
template<typename T>
inline Rcpp::XPtr<ode::Solver<T>> get_solver(SEXP xp) {
  return Rcpp::XPtr<ode::Solver<T>>(xp);
}

// Generic Solver_reset
template<typename SystemType, typename ActiveSystemType>
inline void Solver_reset_impl(SEXP solver_xp, bool active) {
  if (active) {
    get_solver<ActiveSystemType>(solver_xp)->reset();
  } else {
    get_solver<SystemType>(solver_xp)->reset();
  }
}

// Generic Solver_time
template<typename SystemType, typename ActiveSystemType>
inline double Solver_time_impl(SEXP solver_xp, bool active) {
  if (active) {
    return get_solver<ActiveSystemType>(solver_xp)->time();
  } else {
    return get_solver<SystemType>(solver_xp)->time();
  }
}

// Generic Solver_state
template<typename SystemType, typename ActiveSystemType>
inline Rcpp::NumericVector Solver_state_impl(SEXP solver_xp, bool active) {
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
inline Rcpp::NumericVector Solver_times_impl(SEXP solver_xp, bool active) {
  if (active) {
    return Rcpp::wrap(get_solver<ActiveSystemType>(solver_xp)->times());
  } else {
    return Rcpp::wrap(get_solver<SystemType>(solver_xp)->times());
  }
}

// Generic Solver_set_state
template<typename SystemType, typename ActiveSystemType>
inline void Solver_set_state_impl(SEXP solver_xp, Rcpp::NumericVector y, double time, bool active) {
  std::vector<double> yy(y.begin(), y.end());
  if (active) {
    get_solver<ActiveSystemType>(solver_xp)->set_state(yy, time);
  } else {
    get_solver<SystemType>(solver_xp)->set_state(yy, time);
  }
}

// Generic Solver_advance_adaptive
template<typename SystemType, typename ActiveSystemType>
inline void Solver_advance_adaptive_impl(SEXP solver_xp, Rcpp::NumericVector times, bool active) {
  if (active) {
    Rcpp::stop("advance_adaptive() not supported for AD solvers. Use advance_fixed() with pre-computed schedule.");
  }
  std::vector<double> ts(times.begin(), times.end());
  get_solver<SystemType>(solver_xp)->advance_adaptive(ts);
}

// Generic Solver_advance_fixed
template<typename SystemType, typename ActiveSystemType>
inline void Solver_advance_fixed_impl(SEXP solver_xp, Rcpp::NumericVector times, bool active) {
  std::vector<double> ts(times.begin(), times.end());
  if (active) {
    get_solver<ActiveSystemType>(solver_xp)->advance_fixed(ts);
  } else {
    get_solver<SystemType>(solver_xp)->advance_fixed(ts);
  }
}

// Generic Solver_step
template<typename SystemType, typename ActiveSystemType>
inline void Solver_step_impl(SEXP solver_xp, bool active) {
  if (active) {
    Rcpp::stop("step() not supported for AD solvers. Use advance_fixed() instead.");
  }
  get_solver<SystemType>(solver_xp)->step();
}

// Generic Solver_get_collect
template<typename SystemType, typename ActiveSystemType>
inline bool Solver_get_collect_impl(SEXP solver_xp, bool active) {
  if (active) {
    return get_solver<ActiveSystemType>(solver_xp)->get_collect();
  } else {
    return get_solver<SystemType>(solver_xp)->get_collect();
  }
}

// Generic Solver_set_collect
template<typename SystemType, typename ActiveSystemType>
inline void Solver_set_collect_impl(SEXP solver_xp, bool x, bool active) {
  if (active) {
    get_solver<ActiveSystemType>(solver_xp)->set_collect(x);
  } else {
    get_solver<SystemType>(solver_xp)->set_collect(x);
  }
}

// Generic Solver_get_history_size
template<typename SystemType, typename ActiveSystemType>
inline std::size_t Solver_get_history_size_impl(SEXP solver_xp, bool active) {
  if (active) {
    return get_solver<ActiveSystemType>(solver_xp)->get_history_size();
  } else {
    return get_solver<SystemType>(solver_xp)->get_history_size();
  }
}

// Generic Solver_get_history_step
template<typename SystemType, typename ActiveSystemType>
inline Rcpp::DataFrame Solver_get_history_step_impl(SEXP solver_xp, std::size_t i, 
                                                     Rcpp::CharacterVector names, bool active) {
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
inline Rcpp::List Solver_get_history_impl(SEXP solver_xp, Rcpp::CharacterVector names, bool active) {
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
    out[j] = Rcpp::NumericVector(cols[j].begin(), cols[j].end());
  }
  out.attr("names") = names;
  
  return Rcpp::DataFrame(out);
}

// Generic Solver_set_target
template<typename SystemType, typename ActiveSystemType>
inline void Solver_set_target_impl(SEXP solver_xp, 
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
inline Rcpp::List Solver_fit_impl(SEXP solver_xp,
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

#endif // ODELIA_SOLVER_INTERFACE_HPP_