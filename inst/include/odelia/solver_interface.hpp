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
#include <odelia/gradient.hpp>
#include <odelia/rcpp_interface_helpers.hpp>

namespace odelia {
namespace solver {

// Helper to get solver pointer (templated)
template<typename T>
inline Rcpp::XPtr<ode::Solver<T>> get_solver(SEXP xp) {
  return Rcpp::XPtr<ode::Solver<T>>(xp);
}

// The step-wise Solver operations. R holds only the double Solver, so these
// forward straight to it; the AD replay is built internally by the gradient
// drivers at the bottom of this file.

template<typename SystemType>
inline void Solver_reset_impl(SEXP solver_xp) {
  get_solver<SystemType>(solver_xp)->reset();
}

template<typename SystemType>
inline double Solver_time_impl(SEXP solver_xp) {
  return get_solver<SystemType>(solver_xp)->time();
}

template<typename SystemType>
inline Rcpp::NumericVector Solver_state_impl(SEXP solver_xp) {
  return Rcpp::wrap(get_solver<SystemType>(solver_xp)->state());
}

template<typename SystemType>
inline Rcpp::NumericVector Solver_times_impl(SEXP solver_xp) {
  return Rcpp::wrap(get_solver<SystemType>(solver_xp)->times());
}

template<typename SystemType>
inline void Solver_set_state_impl(SEXP solver_xp, Rcpp::NumericVector y, double time) {
  std::vector<double> yy(y.begin(), y.end());
  get_solver<SystemType>(solver_xp)->set_state(yy, time);
}

template<typename SystemType>
inline void Solver_advance_adaptive_impl(SEXP solver_xp, Rcpp::NumericVector times) {
  std::vector<double> ts(times.begin(), times.end());
  get_solver<SystemType>(solver_xp)->advance_adaptive(ts);
}

template<typename SystemType>
inline void Solver_advance_fixed_impl(SEXP solver_xp, Rcpp::NumericVector times) {
  std::vector<double> ts(times.begin(), times.end());
  get_solver<SystemType>(solver_xp)->advance_fixed(ts);
}

template<typename SystemType>
inline void Solver_advance_euler_impl(SEXP solver_xp, Rcpp::NumericVector times) {
  std::vector<double> ts(times.begin(), times.end());
  get_solver<SystemType>(solver_xp)->advance_euler(ts);
}

template<typename SystemType>
inline void Solver_step_impl(SEXP solver_xp) {
  get_solver<SystemType>(solver_xp)->step();
}

template<typename SystemType>
inline bool Solver_get_collect_impl(SEXP solver_xp) {
  return get_solver<SystemType>(solver_xp)->get_collect();
}

template<typename SystemType>
inline void Solver_set_collect_impl(SEXP solver_xp, bool x) {
  get_solver<SystemType>(solver_xp)->set_collect(x);
}

template<typename SystemType>
inline std::size_t Solver_get_history_size_impl(SEXP solver_xp) {
  return get_solver<SystemType>(solver_xp)->get_history_size();
}

template<typename SystemType>
inline Rcpp::DataFrame Solver_get_history_step_impl(SEXP solver_xp, std::size_t i) {
  auto solver = get_solver<SystemType>(solver_xp);
  if (i >= solver->get_history_size()) {
    Rcpp::stop("Index out of bounds");
  }
  Rcpp::CharacterVector names = Rcpp::wrap(solver->get_system().record_colnames());
  std::vector<double> out = solver->get_history_step(i).record_step();

  Rcpp::List df_list(names.size());
  for (size_t j = 0; j < static_cast<size_t>(names.size()); ++j) {
    df_list[j] = out[j];
  }
  df_list.attr("names") = names;
  return Rcpp::DataFrame(df_list);
}

template<typename SystemType>
inline Rcpp::List Solver_get_history_impl(SEXP solver_xp) {
  auto solver = get_solver<SystemType>(solver_xp);
  Rcpp::CharacterVector names = Rcpp::wrap(solver->get_system().record_colnames());
  const int ncols = names.size();
  const size_t nrows = solver->get_history_size();
  std::vector<std::vector<double>> cols(ncols);
  for (auto& col : cols) col.reserve(nrows);

  for (size_t i = 0; i < nrows; ++i) {
    auto row = solver->get_history_step(i).record_step();
    for (int j = 0; j < ncols; ++j) {
      cols[j].push_back(row[j]);
    }
  }

  Rcpp::List out(ncols);
  for (int j = 0; j < ncols; ++j) {
    out[j] = Rcpp::NumericVector(cols[j].begin(), cols[j].end());
  }
  out.attr("names") = names;
  return Rcpp::DataFrame(out);
}

// Marshal an R (observation matrix, 1-based obs_indices) pair into a least_squares
// functional. The functional owns the fit data; the sampling grid is the
// recording's, so no schedule is passed here.
inline ode::least_squares least_squares_from_r(Rcpp::NumericMatrix observations,
                                               Rcpp::IntegerVector obs_indices) {
  ode::least_squares f;

  int nrows = observations.nrow(), ncols = observations.ncol();
  f.observations.resize(nrows);
  for (int i = 0; i < nrows; ++i) {
    f.observations[i].resize(ncols);
    for (int j = 0; j < ncols; ++j) {
      f.observations[i][j] = observations(i, j);
    }
  }

  f.obs_indices.resize(obs_indices.size());
  for (size_t i = 0; i < static_cast<size_t>(obs_indices.size()); ++i) {
    f.obs_indices[i] = obs_indices[i] - 1;  // R 1-based -> C++ 0-based
  }
  return f;
}

// ---- Gradients on the double handle ----------------------------------------
// R holds only the double Solver; these helpers differentiate on the active solver
// (the double System lifted via rebind_from) and return doubles.
template <class SystemType, class ActiveSystemType>
ode::Solver<ActiveSystemType>& active_solver(ode::Solver<SystemType>& d) {
  if (!d.active_solver) {
    d.active_solver = std::make_shared<ode::Solver<ActiveSystemType>>(
        d.get_system_ref().template rebind_from<typename ActiveSystemType::value_type>(),
        d.get_control());
  }
  return *d.active_solver;
}

template <class SystemType, class ActiveSystemType, class Functional>
std::pair<double, std::vector<double>>
gradient_on_double(SEXP solver_xp, const ode::DifferentiationTargets& ind,
                   const std::vector<double>& schedule, Functional&& functional) {
  auto d = get_solver<SystemType>(solver_xp);
  auto& active = active_solver<SystemType, ActiveSystemType>(*d);
  active.set_schedule(schedule);  // hand the recorded L1 schedule to the active twin
  return ode::compute_gradient(active, ind, std::forward<Functional>(functional));
}

template <class SystemType, class ActiveSystemType, class Functional>
std::pair<std::vector<double>, std::vector<std::vector<double>>>
jacobian_on_double(SEXP solver_xp, const ode::DifferentiationTargets& ind,
                   const std::vector<double>& schedule, Functional&& functional) {
  auto d = get_solver<SystemType>(solver_xp);
  auto& active = active_solver<SystemType, ActiveSystemType>(*d);
  active.set_schedule(schedule);  // hand the recorded L1 schedule to the active twin
  return ode::compute_jacobian(active, ind, std::forward<Functional>(functional));
}

// Value + least-squares gradient from one recording, so an optimiser's fn/gr share
// the tape. This is the calibration entry; an arbitrary functional is differentiated
// through gradient_on_double / jacobian_on_double instead (see Solver_gradient_final_state).
template <class SystemType, class ActiveSystemType>
Rcpp::List Solver_value_and_gradient_impl(SEXP solver_xp,
                                          Rcpp::Nullable<Rcpp::NumericVector> ic,
                                          Rcpp::Nullable<Rcpp::NumericVector> params,
                                          Rcpp::NumericVector times,
                                          Rcpp::NumericMatrix observations,
                                          Rcpp::IntegerVector obs_indices) {
  // Seed every param then every IC; `values` follows the same params-then-ics order.
  auto d = get_solver<SystemType>(solver_xp);
  ode::DifferentiationTargets ind;
  if (!params.isNull()) {
    Rcpp::NumericVector v(params);
    for (int i = 0; i < v.size(); ++i) {
      ind.params.push_back(i);
      ind.values.push_back(v[i]);
    }
  }
  if (!ic.isNull()) {
    Rcpp::NumericVector v(ic);
    for (int j = 0; j < v.size(); ++j) {
      ind.ics.push_back(j);
      ind.values.push_back(v[j]);
    }
  }
  // The fit observations are handed in per call and owned by the functional; the
  // solver holds no fit state.
  auto& active = active_solver<SystemType, ActiveSystemType>(*d);
  auto functional = least_squares_from_r(observations, obs_indices);
  // `times` from R is the recorded schedule the active twin replays; least_squares
  // samples the collected trajectory at obs_indices.
  active.set_schedule(std::vector<double>(times.begin(), times.end()));
  auto [value, gradient] = ode::compute_gradient(active, ind, functional);
  return Rcpp::List::create(Rcpp::Named("value") = value,
                            Rcpp::Named("gradient") = Rcpp::wrap(gradient));
}

} // namespace solver
} // namespace odelia

#endif
