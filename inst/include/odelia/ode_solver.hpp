#ifndef ODELIA_ODE_SOLVER_HPP_
#define ODELIA_ODE_SOLVER_HPP_

#include <odelia/ode_solver_internal.hpp>

namespace odelia {
namespace ode {

// This is a wrapper class that is meant to simplify the
// difficuly of ownership semantics around the solver and system.
// It is mostly just be a generic wrapper around ode::Solver<System>
// You would only write your own if you had system-specific needs, e.g. events, cohort introudctions etc.

// TODO
// - add ability to set time directly
// - collect gathers vector of variables at each step
// - move collect into ode::Solver so that it can be used more generally
//   for any system, making this a generic Solver class

template <typename System>
class Solver
{
public:
  Solver(System sys_, OdeControl control) : system(sys_), solver(system, control)
  {
    collect = true;
  }

  // TODO: solver.reset() will set time within the solver to zero.
  // However, there is no other current way of setting the time within
  // the solver.  It might be better to add a set_time method within
  // ode::Solver, and then here do explicitly ode_solver.set_time(0)?
  void reset()
  {
    system.reset();
    solver.reset(system);
    history.clear();
  }

  // collectors
  double time() const { return solver.get_time(); }

  ode::state_type<System> state() const { return solver.get_state(); }
  std::vector<double> times() const { return solver.get_times(); }

  System get_system() const { return system; }
  System& get_system_ref() { return system; }

  void set_state(std::vector<double> y, double time)
  {
    util::check_length(y.size(), system.ode_size());
    internal::set_ode_state(system, y, time);
    solver.reset(system);
    solver.set_state_from_system(system);
  }

  // Take a series of adaptive steps up to some time
  void advance_adaptive(std::vector<double> times)
  {
    if (times.empty())
    {
      util::stop("'times' must be vector of at least length 1");
    }
    std::vector<double>::const_iterator t = times.begin();
    if (!util::identical(*t++, time()))
    {
      util::stop("First element in 'times' must be same as current time");
    }

    if (collect)
    {
      history.push_back(system);
    }

    while (t != times.end())
    {
      solver.advance_adaptive(system, *t++);
      if (collect)
      {
        history.push_back(system);
      }
    }
  }

  // Take a series of steps at specified time steps
  void advance_fixed(std::vector<double> times)
  {
    if (times.empty())
    {
      util::stop("'times' must be vector of at least length 1");
    }
    std::vector<double>::const_iterator t = times.begin();
    if (!util::identical(*t++, time()))
    {
      util::stop("First element in 'times' must be same as current time");
    }

    if (collect)
    {
      history.push_back(system);
    }

    while (t != times.end())
    {
      solver.step_to(system, *t++);
      if (collect)
      {
        history.push_back(system);
      }
    }
  }

  void step()
  {
    solver.step(system);
    if (collect)
    {
      history.push_back(system);
    }
  }

  bool get_collect() const { return collect; }

  void set_collect(bool x) { collect = x; }

  std::size_t get_history_size() const { return history.size(); }

  std::vector<System> get_history() const { return history; }

  System get_history_step(std::size_t i) const { return history.at(i); }

  // Fit configuration methods
  void set_target(const std::vector<double>& times, 
                 const std::vector<std::vector<double>>& targets,
                 const std::vector<size_t>& obs_indices) {
    fit_times_ = times;
    targets_ = targets;
    obs_indices_ = obs_indices;
  }
  // Advance solver and return states only at observation times
  std::vector<std::vector<typename System::value_type>> advance_target() {
    // Check that target has been set
    if (fit_times_.empty()) {
      util::stop("Must call set_target() before advance_target()");
    }
    
    // Check starting time matches
    if (!util::identical(fit_times_[0], time())) {
      util::stop("First element in fit_times must be same as current time");
    }
    
    // Vector to store states at observation times
    std::vector<std::vector<typename System::value_type>> observations;
    observations.reserve(obs_indices_.size());
    
    // Track which observation we're looking for
    size_t obs_idx = 0;
    
    // Check if initial time is an observation
    if (obs_idx < obs_indices_.size() && obs_indices_[obs_idx] == 0) {
      observations.push_back(state());
      obs_idx++;
    }
    
    // Step through times
    for (size_t i = 1; i < fit_times_.size(); ++i) {
      solver.step_to(system, fit_times_[i]);
      
      // Check if this time index is an observation point
      while (obs_idx < obs_indices_.size() && obs_indices_[obs_idx] == i) {
        observations.push_back(state());
        obs_idx++;
      }
    }
    
    return observations;
  }
  const std::vector<double>& fit_times() const { return fit_times_; }
  const std::vector<std::vector<double>>& targets() const { return targets_; }
  const std::vector<size_t>& obs_indices() const { return obs_indices_; }

  // Should we record history at every step?
  // TODO: should this be part of ode_solver?
  std::vector<System> history;

private:
  bool collect;
  System system;
  SolverInternal<System> solver;
  
  // Fit configuration for AD gradient computation
  std::vector<double> fit_times_;
  std::vector<size_t> obs_indices_;
  std::vector<std::vector<double>> targets_;

};
}
}
#endif
