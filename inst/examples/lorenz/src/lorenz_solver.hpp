#include <odelia/ode_solver.hpp>
#include "lorenz_system.hpp"

namespace odelia {
namespace ode {

// This is a little wrapper class that is meant to simplify the
// difficuly of ownership semantics around the solver and object.
class Runner {

  public:
  Runner(LorenzSystem sys_, ode::OdeControl control) : sys(sys_), solver(sys, control) {}
  
  double time() const {return solver.get_time();}
  
  ode::state_type state() const {return solver.get_state();}
  
  void set_state(ode::state_type y, double time) {
    util::check_length(y.size(), sys.ode_size());
    ode::internal::set_ode_state(sys, y, time);
    solver.reset(sys);
    solver.set_state_from_system(sys);
  }
  
  void set_state_from_system() {
    solver.set_state_from_system(sys);
  }

  std::vector<double> times() const {
    return solver.get_times();
  }
  
  LorenzSystem system() const {return sys;}

  void advance_adaptive(double time) { 
    solver.advance_adaptive(sys, time); 
  }
  
  void advance_fixed(std::vector<double> times) {
    solver.advance_fixed(sys, times);
  }
  void step() {
    solver.step(sys);
  }
  
  void step_to(double time) {
    solver.step_to(sys, time);
  }

  LorenzSystem sys;
  ode::Solver<LorenzSystem> solver;
};

}
}
