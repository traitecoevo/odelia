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
#include <string>
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

// The parameters of the System this solver steps. A Solver copies the System it
// is built from and steps the copy, so parameters set on the original afterwards
// do not reach it; these read and write the copy's, in the order ad_parameters()
// lists them, which is the order set_params() reads and a gradient comes in.
template<typename SystemType>
inline Rcpp::NumericVector Solver_params_impl(SEXP solver_xp) {
  SystemType& sys = get_solver<SystemType>(solver_xp)->get_system_ref();
  std::vector<double> out;
  for (double* p : sys.ad_parameters()) out.push_back(*p);
  return Rcpp::wrap(out);
}

// Re-seeds the solver at its current state and time, as set_state does, so the
// rates it carries are the new parameters' and the step size starts afresh.
template<typename SystemType>
inline void Solver_set_params_impl(SEXP solver_xp, Rcpp::NumericVector params) {
  auto solver = get_solver<SystemType>(solver_xp);
  SystemType& sys = solver->get_system_ref();
  if (static_cast<std::size_t>(params.size()) != sys.ad_parameters().size()) {
    Rcpp::stop("'params' needs one entry per parameter");
  }
  std::vector<double> pv(params.begin(), params.end());
  sys.set_params(pv.begin());
  solver->set_state(solver->state(), solver->time());
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
// The forward pass steps `schedule` -- the reference run's own step times -- at
// the parameters and initial state given, keeping its states, so the loss is a
// function of those inputs alone and not of a step-size controller that would
// move with them. The loss is sum over observations i and state components j of
// (y_j(schedule[obs_indices[i]]) - target(i, j))^2.
//
// The gradient comes one of two ways, `mode`, over the same steps under the
// solver's own method:
//
//   "reverse"  one sweep over the recording, taken interval by interval between
//              observed rows: each observation's own term is added to the
//              adjoint at its row, and the sweep carries the sum down to the
//              next one. What arrives at row 0 is dL/dy0; what the sweep
//              accumulates is dL/dtheta. Its cost does not grow with the count
//              of inputs. RODAS cannot be swept.
//   "forward"  the schedule walked once more at a scalar carrying a direction
//              per input, W at a time (sweep_tangent), each observation's term
//              taken against dy/dtheta at its row. Its cost grows with the count
//              of inputs, and it runs under every method, RODAS included.
//
// An empty `mode` is "forward" under RODAS and "reverse" otherwise.
//
// Returns list(loss, gradient) with the gradient's entries in the order the old
// fitting interface used: the parameters (as the System's ad_parameters() lists
// them) if `params` was given, then the initial state if `initial_state` was.
template<typename SystemType, std::size_t W = 8>
inline Rcpp::List Solver_loss_and_gradient_impl(
    SEXP solver_xp, Rcpp::NumericVector schedule, Rcpp::NumericMatrix target,
    Rcpp::IntegerVector obs_indices,
    Rcpp::Nullable<Rcpp::NumericVector> initial_state,
    Rcpp::Nullable<Rcpp::NumericVector> params, std::string mode = "") {
  if (initial_state.isNull() && params.isNull()) {
    Rcpp::stop("Must provide at least one of 'initial_state' or 'params'");
  }
  if (obs_indices.size() == 0) {
    Rcpp::stop("'obs_indices' must name at least one observation");
  }
  const std::size_t n_times = schedule.size();
  if (n_times < 2) {
    Rcpp::stop("'schedule' must hold the run's start and at least one step");
  }
  for (std::size_t k = 1; k < n_times; ++k) {
    if (!(schedule[k] > schedule[k - 1])) {
      Rcpp::stop("'schedule' must be strictly increasing");
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
  if (initial_state.isNotNull()) {
    Rcpp::NumericVector y(initial_state);
    if (static_cast<std::size_t>(y.size()) != n) {
      Rcpp::stop("'initial_state' needs one entry per state variable");
    }
    std::vector<double> yv(y.begin(), y.end());
    sys.set_initial_state(yv.begin(), schedule[0]);
  }
  sys.reset();
  std::vector<double> y0(n);
  sys.ode_state(y0.begin());

  std::vector<ode::schedule_row> rows(n_times);
  rows[0] = {schedule[0], std::numeric_limits<double>::quiet_NaN()};
  for (std::size_t k = 1; k < n_times; ++k) {
    rows[k] = {schedule[k], schedule[k] - schedule[k - 1]};
  }

  // Whether this System can be walked forward at all, and under RODAS: known at
  // compile time, so a System that cannot is never instantiated at the lane
  // scalar, and keeps the reverse sweep it had.
  using lane_system =
      typename ode::rebound_system<SystemType, ode::tangent_lanes<W>>::type;
  constexpr bool walkable = ode::TangentSweepable<SystemType, ode::tangent_lanes<W>>;
  const ode::Method method = solver->method();
  if (mode.empty()) {
    mode = method == ode::Method::rodas ? "forward" : "reverse";
  }
  if (mode != "forward" && mode != "reverse") {
    Rcpp::stop("'mode' must be \"forward\" or \"reverse\"");
  }
  if (method == ode::Method::ros34pw2) {
    Rcpp::stop("method = 'ros34pw2' cannot be differentiated: the steps at "
               "which it formed a new Jacobian are not recorded, so no walk "
               "of the run takes them again; use a solver built with 'rodas', "
               "'rkck' or 'dopri'");
  }
  if (mode == "reverse" && method == ode::Method::rodas) {
    Rcpp::stop("method = 'rodas' cannot be swept in reverse; use mode = "
               "\"forward\", or a solver built with 'rkck' or 'dopri'");
  }
  std::vector<std::size_t> observed(obs_indices.size());
  for (int i = 0; i < obs_indices.size(); ++i) {
    const int idx = obs_indices[i];
    if (idx < 1 || static_cast<std::size_t>(idx) > n_times) {
      Rcpp::stop("'obs_indices' must index into 'schedule' (1-based)");
    }
    observed[i] = static_cast<std::size_t>(idx - 1);
  }

  if (mode == "forward") {
    if constexpr (!walkable) {
      Rcpp::stop("this System cannot be differentiated in forward mode (see "
                 "TangentSweepable)%s",
                 method == ode::Method::rodas
                     ? ", and a run under method = 'rodas' cannot be swept in "
                       "reverse; use a solver built with 'rkck' or 'dopri'"
                     : "; use mode = \"reverse\"");
    } else {
    if (method == ode::Method::rodas &&
        !ode::RodasStep<lane_system>::supported) {
      Rcpp::stop("under method = 'rodas' this System is differentiated in "
                 "forward mode, which needs it to declare `static constexpr "
                 "bool nests_tangents = true` (see NestsTangents) or have its "
                 "own ode_jacobian(); or use a solver built with 'rkck' or "
                 "'dopri'");
    }
    ode::tangent_directions wrt;
    if (params.isNotNull()) {
      for (std::size_t q = 0; q < n_par; ++q) wrt.parameters.push_back(q);
    }
    if (initial_state.isNotNull()) {
      for (std::size_t j = 0; j < n; ++j) wrt.initial.push_back(j);
    }
    const std::vector<ode::tangent_row> d = ode::sweep_tangent<W>(
        sys, solver->control(), method, rows, wrt, observed);
    double loss = 0.0;
    std::vector<double> gradient(wrt.size(), 0.0);
    for (std::size_t i = 0; i < observed.size(); ++i) {
      for (std::size_t j = 0; j < n; ++j) {
        const double diff = d[i].value[j] - target(i, j);
        loss += diff * diff;
        for (std::size_t k = 0; k < gradient.size(); ++k) {
          gradient[k] += 2.0 * diff * d[i].at(j, k);
        }
      }
    }
    return Rcpp::List::create(Rcpp::Named("loss") = loss,
                              Rcpp::Named("gradient") = Rcpp::wrap(gradient));
    }
  }

  ode::Solver<SystemType> replay(sys, solver->control(), method);
  replay.set_collect(false);
  replay.set_keep_states(true);
  replay.set_state(y0, schedule[0]);
  replay.advance_schedule(rows);
  const auto rec = replay.recording();
  // One row per scheduled step, or a step was subdivided and its row no longer
  // describes one Runge-Kutta step -- which the sweep would transpose as one.
  if (rec.size() != n_times) {
    Rcpp::stop("the run did not take the schedule's steps one row each, so it "
               "cannot be swept");
  }

  // Each observed row's term, summed where two observations share a row.
  std::map<std::size_t, ode::adjoint_batch> adjoints;
  double loss = 0.0;
  for (std::size_t i = 0; i < observed.size(); ++i) {
    const std::size_t row = observed[i];
    auto found = adjoints.find(row);
    if (found == adjoints.end()) {
      found = adjoints.emplace(row, ode::adjoint_batch(1, n)).first;
    }
    for (std::size_t j = 0; j < n; ++j) {
      const double diff = rec[row].state[j] - target(i, j);
      loss += diff * diff;
      found->second[0][j] += 2.0 * diff;
    }
  }
  std::vector<std::size_t> at_rows;
  std::vector<ode::adjoint_batch> at;
  for (const auto& [row, adjoint] : adjoints) {
    at_rows.push_back(row);
    at.push_back(adjoint);
  }
  ode::adjoint_batch lambda;
  ode::adjoint_batch parameter_adjoint(1, n_par);
  replay.sweep_adjoint_at(at_rows, at, lambda, parameter_adjoint);

  std::vector<double> gradient;
  if (params.isNotNull()) {
    for (std::size_t q = 0; q < n_par; ++q) gradient.push_back(parameter_adjoint[0][q]);
  }
  if (initial_state.isNotNull()) {
    for (std::size_t j = 0; j < n; ++j) gradient.push_back(lambda[0][j]);
  }
  return Rcpp::List::create(Rcpp::Named("loss") = loss,
                            Rcpp::Named("gradient") = Rcpp::wrap(gradient));
}


// The state at each of `times` and its derivative with respect to the chosen
// parameters (`parameters`, 0-based into ad_parameters()) and components of the
// state at times[0] (`initial`, 0-based), by forward mode.
//
// The solve starts from where this solver stands, with its step size afresh --
// a copy of it is advanced, and the solver itself does not move -- and lands a
// step on every time. Its recording is then walked again
// carrying the directions (sweep_tangent), under the solver's own method.
//
// Returns list(time, y, dy): y is n_times x n, dy an n_times x n x P array, the
// parameters' directions first, then the initial state's, each in the order
// given.
template<typename SystemType, std::size_t W = 8>
inline Rcpp::List Solver_sensitivities_impl(SEXP solver_xp,
                                            Rcpp::NumericVector times,
                                            Rcpp::IntegerVector parameters,
                                            Rcpp::IntegerVector initial) {
  auto solver = get_solver<SystemType>(solver_xp);
  if (solver->method() == ode::Method::ros34pw2) {
    Rcpp::stop("method = 'ros34pw2' cannot be differentiated: the steps at "
               "which it formed a new Jacobian are not recorded, so no walk "
               "of the run takes them again; use a solver built with 'rodas', "
               "'rkck' or 'dopri'");
  }
  const std::size_t n_t = times.size();
  if (n_t == 0) {
    Rcpp::stop("'times' must hold at least the current time");
  }
  // A copy that has stepped keeps no record of where it stands; set_state
  // starts a new one there (and the step size afresh).
  ode::Solver<SystemType> run(*solver);
  run.set_collect(false);
  run.set_state(run.state(), run.time());
  run.set_keep_states(true);
  std::vector<std::size_t> rows{0};
  for (std::size_t k = 1; k < n_t; ++k) {
    run.advance_adaptive(std::vector<double>{times[k - 1], times[k]});
    rows.push_back(run.recording().size() - 1);
  }
  if (n_t == 1 && !util::identical(times[0], run.time())) {
    Rcpp::stop("First element in 'times' must be same as current time");
  }

  ode::tangent_directions wrt;
  for (const int p : parameters) {
    if (p < 0) Rcpp::stop("'parameters' must be non-negative");
    wrt.parameters.push_back(static_cast<std::size_t>(p));
  }
  for (const int i : initial) {
    if (i < 0) Rcpp::stop("'initial' must be non-negative");
    wrt.initial.push_back(static_cast<std::size_t>(i));
  }
  const std::vector<ode::tangent_row> d =
      run.template sweep_tangent<W>(wrt, rows);

  const std::size_t n = d.front().value.size();
  const std::size_t n_dir = wrt.size();
  Rcpp::NumericMatrix y(n_t, n);
  Rcpp::NumericVector dy(n_t * n * n_dir);
  for (std::size_t k = 0; k < n_t; ++k) {
    if (d[k].value.size() != n) {
      Rcpp::stop("the state changed width between the first and last time");
    }
    for (std::size_t i = 0; i < n; ++i) {
      y(k, i) = d[k].value[i];
      for (std::size_t j = 0; j < n_dir; ++j) {
        dy[k + n_t * (i + n * j)] = d[k].at(i, j);
      }
    }
  }
  dy.attr("dim") = Rcpp::IntegerVector::create(n_t, n, n_dir);
  return Rcpp::List::create(Rcpp::Named("time") = times,
                            Rcpp::Named("y") = y,
                            Rcpp::Named("dy") = dy);
}

} // namespace solver
} // namespace odelia

#endif
