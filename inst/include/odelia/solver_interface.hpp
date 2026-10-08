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
#include <limits>
#include <map>
#include <vector>
#include <XAD/XAD.hpp>
#include <odelia/ode_solver.hpp>
#include <odelia/rcpp_interface_helpers.hpp>

namespace odelia {
namespace solver {

// Helper to get solver pointer (templated)
template<typename T>
inline Rcpp::XPtr<ode::Solver<T>> get_solver(SEXP xp) {
  return Rcpp::XPtr<ode::Solver<T>>(xp);
}

// The step-wise Solver operations. R holds only the double Solver, so these
// forward straight to it.

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

// A least-squares loss against observations of the state, and its exact gradient
// with respect to the System's parameters and/or its initial state.
//
// The forward pass replays the schedule `times` records -- the reference run's own
// step times -- at the parameters and initial state given, keeping its states, so
// the loss is a function of those inputs alone and not of a step-size controller
// that would move with them. The loss is sum over observations i and state
// components j of (y_j(times[obs_indices[i]]) - target(i, j))^2.
//
// The gradient is one reverse sweep over that recording, taken as one range per
// interval between observed rows: each observation's own term is added to the
// adjoint at its row, and the sweep carries the sum down to the next one. What
// arrives at row 0 is dL/dy0; what the sweep accumulates is dL/dtheta.
//
// Returns list(loss, gradient) with the gradient's entries in the order the old
// fitting interface used: the parameters (as the System's ad_parameters() lists
// them) if `params` was given, then the initial state if `ic` was.
template<typename SystemType>
inline Rcpp::List Solver_fit_impl(SEXP solver_xp,
                                  Rcpp::NumericVector times,
                                  Rcpp::NumericMatrix target,
                                  Rcpp::IntegerVector obs_indices,
                                  Rcpp::Nullable<Rcpp::NumericVector> ic,
                                  Rcpp::Nullable<Rcpp::NumericVector> params) {
  if (ic.isNull() && params.isNull()) {
    Rcpp::stop("Must provide at least one of 'ic' or 'params'");
  }
  const std::size_t n_times = times.size();
  if (n_times < 2) {
    Rcpp::stop("'times' must hold the run's start and at least one step");
  }
  for (std::size_t k = 1; k < n_times; ++k) {
    if (!(times[k] > times[k - 1])) {
      Rcpp::stop("'times' must be strictly increasing");
    }
  }
  if (target.nrow() != obs_indices.size()) {
    Rcpp::stop("'target' needs one row per entry of 'obs_indices'");
  }

  auto solver = get_solver<SystemType>(solver_xp);
  SystemType sys = solver->get_system();
  const std::size_t n = sys.ode_size();
  const std::size_t n_par = sys.ad_parameters().size();
  if (static_cast<std::size_t>(target.ncol()) != n) {
    Rcpp::stop("'target' needs one column per state variable");
  }
  if (params.isNotNull()) {
    Rcpp::NumericVector p(params);
    if (static_cast<std::size_t>(p.size()) != n_par) {
      Rcpp::stop("'params' needs one entry per parameter");
    }
    std::vector<double> pv(p.begin(), p.end());
    sys.set_params(pv.begin());
  }
  if (ic.isNotNull()) {
    Rcpp::NumericVector y(ic);
    if (static_cast<std::size_t>(y.size()) != n) {
      Rcpp::stop("'ic' needs one entry per state variable");
    }
    std::vector<double> yv(y.begin(), y.end());
    sys.set_initial_state(yv.begin(), times[0]);
  }
  sys.reset();
  std::vector<double> y0(n);
  sys.ode_state(y0.begin());

  std::vector<ode::instruction> program(n_times);
  program[0] = {times[0], std::numeric_limits<double>::quiet_NaN()};
  for (std::size_t k = 1; k < n_times; ++k) {
    program[k] = {times[k], times[k] - times[k - 1]};
  }

  ode::Solver<SystemType> replay(sys, solver->control());
  replay.set_collect(false);
  replay.set_keep_states(true);
  replay.set_state(y0, times[0]);
  replay.advance_recorded(program);
  const auto rec = replay.recording();
  // One row per scheduled step, or a step was subdivided and its row no longer
  // describes one Runge-Kutta step -- which the sweep would transpose as one.
  if (rec.size() != n_times) {
    Rcpp::stop("the replay did not take the schedule's steps one row each, so it "
               "cannot be swept");
  }

  // Each observed row's term, summed where two observations share a row.
  std::map<std::size_t, std::vector<double>> seeds;
  double loss = 0.0;
  for (int i = 0; i < obs_indices.size(); ++i) {
    const int idx = obs_indices[i];
    if (idx < 1 || static_cast<std::size_t>(idx) > n_times) {
      Rcpp::stop("'obs_indices' must index into 'times' (1-based)");
    }
    const std::size_t row = static_cast<std::size_t>(idx - 1);
    std::vector<double>& seed = seeds[row];
    seed.resize(n, 0.0);
    for (std::size_t j = 0; j < n; ++j) {
      const double diff = rec[row].state[j] - target(i, j);
      loss += diff * diff;
      seed[j] += 2.0 * diff;
    }
  }

  ode::adjoint_rows lambda(1, n);
  ode::adjoint_rows parameter_adjoint(1, n_par);
  auto it = seeds.rbegin();
  std::size_t hi = it->first;
  for (std::size_t j = 0; j < n; ++j) lambda[0][j] = it->second[j];
  for (++it; it != seeds.rend(); ++it) {
    replay.solve_adjoint(lambda, parameter_adjoint, it->first, hi);
    for (std::size_t j = 0; j < n; ++j) lambda[0][j] += it->second[j];
    hi = it->first;
  }
  if (hi > 0) {
    replay.solve_adjoint(lambda, parameter_adjoint, 0, hi);
  }

  std::vector<double> gradient;
  if (params.isNotNull()) {
    for (std::size_t q = 0; q < n_par; ++q) gradient.push_back(parameter_adjoint[0][q]);
  }
  if (ic.isNotNull()) {
    for (std::size_t j = 0; j < n; ++j) gradient.push_back(lambda[0][j]);
  }
  return Rcpp::List::create(Rcpp::Named("loss") = loss,
                            Rcpp::Named("gradient") = Rcpp::wrap(gradient));
}

} // namespace solver
} // namespace odelia

#endif
