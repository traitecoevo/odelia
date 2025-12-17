#include <odelia/ode_solver.hpp>
#include "lorenz_system.hpp"

namespace odelia {
namespace ode {

// This is a wrapper class that is meant to simplify the
// difficuly of ownership semantics around the solver and system.
// It is mostly just be a generic wrapper around ode::Solver<System>

// TODO
// - use templates to allow this to be generic for any System
// - move collect into ode::Solver so that it can be used more generally
//   for any system, making this a generic Runner class 
// You would only write your own if you had system-specific needs, e.g. events, cohort introudctions etc.

class Runner {

  public:

  Runner(LorenzSystem sys_, ode::OdeControl control) : sys(sys_), solver(sys, control) {   
    collect = true;
  }

  // TODO: solver.reset() will set time within the solver to zero.
  // However, there is no other current way of setting the time within
  // the solver.  It might be better to add a set_time method within
  // ode::Solver, and then here do explicitly ode_solver.set_time(0)?
  void reset()
  {
    sys.reset();
    solver.reset(sys);
    history.clear();
  }

  // collectors
  double time() const {return solver.get_time();}
  
  ode::state_type state() const {return solver.get_state();}
  std::vector<double> times() const { return solver.get_times(); }

  LorenzSystem system() const { return sys; }

  void set_state(ode::state_type y, double time) {
    util::check_length(y.size(), sys.ode_size());
    ode::internal::set_ode_state(sys, y, time);
    solver.reset(sys);
    solver.set_state_from_system(sys);
  }

  // Take a series of adaptive steps up to some time
  void advance_adaptive(std::vector<double> times) {
    if (times.empty()) {
      util::stop("'times' must be vector of at least length 1");
    }
    std::vector<double>::const_iterator t = times.begin();
    if (!util::identical(*t++, time())) {
      util::stop("First element in 'times' must be same as current time");
    }
    
    if (collect) { history.push_back(sys); }
    
    while (t != times.end()) {
      solver.advance_adaptive(sys, *t++);
      if (collect) { history.push_back(sys); }
    }
  }

  // Take a series of steps at specified time steps
  void advance_fixed(std::vector<double> times) {
    if (times.empty()) {
      util::stop("'times' must be vector of at least length 1");
    }
    std::vector<double>::const_iterator t = times.begin();
    if (!util::identical(*t++, time())) {
      util::stop("First element in 'times' must be same as current time");
    }
    
    if (collect) { history.push_back(sys); }
    
    while (t != times.end()) {
      solver.step_to(sys, *t++);
      if (collect) { history.push_back(sys); }
    }
  }

  void step() {
    solver.step(sys);
    if (collect) { history.push_back(sys); }
  } 

  
  bool get_collect() const { return collect; }

  void set_collect(bool x) { collect = x; }

  std::vector<LorenzSystem> get_history() const { return history; }

  LorenzSystem get_history_element(std::size_t i) const { return history.at(i);}
  
  std::size_t get_history_size() const { return history.size(); }

private:
  // Should we record history at every step?
  // TODO: should this be part of ode_solver?
  bool collect;
  std::vector<LorenzSystem> history;

  LorenzSystem sys;
  ode::Solver<LorenzSystem> solver;
};

}
}
