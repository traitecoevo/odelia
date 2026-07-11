// -*-c++-*-
#ifndef ODELIA_ODE_SOLVER_INTERNAL_HPP_
#define ODELIA_ODE_SOLVER_INTERNAL_HPP_

#include <odelia/ode_interface.hpp>
#include <odelia/ode_control.hpp>
#include <odelia/ode_step.hpp>
#include <odelia/ode_step_rodas.hpp>

#include <limits>
#include <vector>
#include <cstddef>

namespace odelia {
namespace ode {

// Integration method: the explicit Cash-Karp RKCK 4(5) stepper (default) or the
// implicit RODAS4(3) Rosenbrock stepper for stiff systems.
enum class Method { rkck, rodas };

template <class System>
class SolverInternal {
public:
  // Extract scalar type from System using traits
  using value_type = typename System::value_type;
  using state_type = std::vector<value_type>;

  SolverInternal(const System &system, OdeControl control_,
                 Method method_ = Method::rkck);
  void reset(const System& system);
  void set_state_from_system(const System& system);

  state_type get_state() const {return y;}
  double get_time() const {return time;}
  std::vector<double> get_times() const {return prev_times;}

  void advance_adaptive(System &system, double time_max_);
  void advance_fixed(System& system, const std::vector<double>& times);
  void advance_euler(System& system, const std::vector<double>& times);

  void step(System& system);
  void step_to(System& system, double time_max_);
  void step_euler(System& system, double time_max_);

  void set_time_max(double time_max_);

private:
  void resize(size_t size_);
  void setup_dydt_in(System& system);
  void save_dydt_out_as_in();
  void set_time(double t);

  // Stepper dispatch: SolverInternal holds both steppers and forwards to the one
  // selected at construction. The adaptive controller (see step()) is otherwise
  // stepper-agnostic.
  void stepper_step(System& system, double time_, double step_size,
                    state_type& y_, state_type& yerr_,
                    const state_type& dydt_in_, state_type& dydt_out_) {
    if (method == Method::rodas) {
      if constexpr (RodasStep<System>::supported) {
        rodas_stepper.step(system, time_, step_size, y_, yerr_, dydt_in_,
                           dydt_out_);
      } else {
        // RODAS is unavailable for this System: either it provides no rebind_from()
        // hook for the AD Jacobian, or its scalar type is itself active (nested
        // tangent-over-adjoint is not yet wired up -- see issue #35). The passive
        // solver of a system with rebind_from() supports RODAS.
        util::stop("method='rodas' is not available for this system/scalar type "
                   "(needs a rebind_from() hook and a non-active scalar); "
                   "use method='rkck'.");
      }
    } else {
      stepper.step(system, time_, step_size, y_, yerr_, dydt_in_, dydt_out_);
    }
  }
  size_t stepper_order() const {
    return method == Method::rodas ? rodas_stepper.order() : stepper.order();
  }
  bool stepper_can_use_dydt_in() const {
    return method == Method::rodas ? RodasStep<System>::can_use_dydt_in
                                   : Step<System>::can_use_dydt_in;
  }
  bool stepper_first_same_as_last() const {
    return method == Method::rodas ? RodasStep<System>::first_same_as_last
                                   : Step<System>::first_same_as_last;
  }

  OdeControl control;
  Method method;
  Step<System> stepper;
  RodasStep<System> rodas_stepper;

  double step_size_last; // Size of last successful step (or suggestion)

  double time;     // Current time
  double time_max; // Time we will not go past
  std::vector<double> prev_times; // Vector of previous times.

  state_type y;        // Vector of current system state
  state_type yerr;     // Vector of error estimates
  state_type dydt_in;  // Vector of dydt at beginning of step
  state_type dydt_out; // Vector of dydt during step

  bool dydt_in_is_clean;
};

// NOTE I'm setting the initial system size to 0 here, but some
// systems are self-initialising.
template <class System>
SolverInternal<System>::SolverInternal(const System &system, OdeControl control_,
                                       Method method_)
  : control(control_), method(method_) {
  reset(system);
}

// NOTE: This resets *everything* to basically a recreated object.
template <class System>
void SolverInternal<System>::reset(const System& system) {
  prev_times.clear();
  step_size_last = control.step_size_initial;
  time_max = std::numeric_limits<double>::infinity();
  set_state_from_system(system);
}

// Record this ODE step's node positions on a Replayable System during the adaptive
// pass; a no-op for any System that doesn't record.
template <typename System>
void record_ode_step(System& system) {
  if constexpr (Replayable<System>) {
    system.record_ode_step();
  }
}

// On the replay pass, let a Replayable System restore what it recorded for this step;
// a no-op otherwise. Called from step_to. What restoring does -- reuse recorded field
// values, or nothing so the field is recomputed with the active scalar -- is the
// System's own choice.
template <typename System>
void replay_step(System& system) {
  if constexpr (Replayable<System>) {
    system.replay_step();
  }
}

template <class System>
void SolverInternal<System>::set_state_from_system(const System& system) {
  set_time(ode::ode_time(system));
  resize(system.ode_size());
  system.ode_state(y.begin());
  system.ode_rates(dydt_in.begin());
  dydt_in_is_clean = true;
}

template <class System>
void SolverInternal<System>::advance_adaptive(System &system, double time_max_)
{
  set_time_max(time_max_);
  while (time < time_max) {
    step(system);
  }
}

// NOTE: We take a vector of times {t_0, t_1, ...}.  This vector
// *must* contain a starting time, but can otherwise be empty.  We
// will step exactly to t_1, then to t_2 up to the end point.  No step
// size adjustments will be done.  This is used in the SCM.
//
// NOTE: Careful here: exact floating point comparison in determining
// that we're starting from the right place.  However, because we take
// care to return and add end points exactly, this should actually be
// the correct move.
template <class System>
void SolverInternal<System>::advance_fixed(System& system,
                                   const std::vector<double>& times) {
  if (times.empty()) {
    util::stop("'times' must be vector of at least length 1");
  }
  std::vector<double>::const_iterator t = times.begin();
  if (!util::identical(*t++, time))
  {
    util::stop("First element in 'times' must be same as current time");
  }
  while (t != times.end()) {
    step_to(system, *t++);
  }
}

// Plain forward (explicit) Euler integration over a supplied grid
// {t_0, t_1, ...}.  Unlike advance_fixed (which still drives the full multi-stage
// RKCK stepper at each interval), this does ONE derivative evaluation per
// interval: derivatives at the current state, then y <- y + h * dydt, advancing
// the time exactly to each grid point.  The `Step` (RKCK) machinery is bypassed
// entirely, so there is no error estimate and no step-size control.  Used to run
// systems the way fixed-step DGVMs do (e.g. a daily step).
template <class System>
void SolverInternal<System>::advance_euler(System& system,
                                   const std::vector<double>& times) {
  if (times.empty()) {
    util::stop("'times' must be vector of at least length 1");
  }
  std::vector<double>::const_iterator t = times.begin();
  if (!util::identical(*t++, time)) {
    util::stop("First element in 'times' must be same as current time");
  }
  while (t != times.end()) {
    step_euler(system, *t++);
  }
}

// A single forward-Euler step from the current time up to time_max_.  One
// derivative evaluation; no error estimate.  Leaves the system synchronised with
// the new state (like step_to, whose final RK derivs settles the system at y) so
// that collected history / record_step reflect the post-step values.
template <class System>
void SolverInternal<System>::step_euler(System& system, double time_max_) {
  set_time_max(time_max_);
  const double h = time_max - time;
  // Derivatives at the current state (also sets the system to y at this time).
  ode::derivs(system, y, dydt_in, time);
  const size_t size = y.size();
  for (size_t i = 0; i < size; ++i) {
    y[i] += h * dydt_in[i];
  }
  time = time_max;
  // Settle the system onto the new state at the new time.
  ode::internal::set_ode_state(system, y, time);
  prev_times.push_back(time);
  dydt_in_is_clean = false;
}

// After `stepper.step()`, the GSL checks to see if the step succeeded
// (some steppers look like they fail for non-user function error),
// and the divides the step size by 2.  If it fails with `EFAULT` or
// `EBADFUNC`, then it aborts.  The only place that errors are
// actually checked in the user function, and the two errors that
// cause abort are the only two that should be thrown there.
//
// There are several different logical step sizes:
//
// 1. this->step_size_last: Size of the last successful step last
//    time, or a suggestion of one.  This will get updated as leave
//    the function only if (1) the step is successful and (2) if we're
//    not in the final step.  It's not actually quite the size of the
//    last step, either -- it's the size that the controller suggested
//    updating the step size too after the last current step.
//
// 2. step_size: The size that the current iteration actually advanced
//    the system (or will) via `stepper.step`.
//
// 3. step_size_next: The size of the proposed next step (or retry of
//    the current step).
template <class System>
void SolverInternal<System>::step(System& system) {
  const double time_orig = time, time_remaining = time_max - time;
  double step_size = step_size_last;


  // Save y in case of failure in a step (recall that stepper.step
  // changes 'y')
  const state_type y_orig = y;
  const size_t size = y.size();

  // Compute the derivatives at the beginning.
  setup_dydt_in(system);

  while (true) {
    // Does this appear to be the last step before reaching `time_max`?
    const bool final_step = step_size > time_remaining;
    if (final_step) {
      step_size = time_remaining;
    }

    stepper_step(system, time, step_size, y, yerr, dydt_in, dydt_out);

    const double step_size_next =
      control.adjust_step_size(size, stepper_order(), step_size,
			       y, yerr, dydt_out);

    if (control.step_size_shrank()) {
        // GSL checks that the step size is actually decreased.
        // Probably we can do this by comparing against hmin?  There are
        // probably loops that this will not catch, but require that
        // hmin << t
         const double time_next = time + step_size_next;
      if (step_size_next < step_size && time_next > time_orig) {
      	// Step was decreased. Undo step (resetting the state y and
        // time), and try again with the new step_size.
      	y         = y_orig;
      	time      = time_orig;
      	step_size = step_size_next;
      } else {
      	// We've reached limits of machine accuracy in differences of
      	// step sizes or time (or both).
        util::stop("Cannot achieve the desired accuracy");
      }
    } else {
      // We have successfully taken a step and will return.  Update
      // time to reflect this, ensuring that if we're on the last step
      // we will end up exactly at time_max.
      //
      // Suggest step size for next time-step. Change of step size is not
      //  suggested in the final step, because that step can be very
      //  small compared to previous step, to reach time_max.
      if (final_step) {
	      time = time_max;
      } else {
	      time += step_size;
	      step_size_last = step_size_next;
      }
      prev_times.push_back(time);
      save_dydt_out_as_in();
      record_ode_step(system);
      return; // This exits the infinite loop.
    }
  }
}

// This takes a step up to time "time_max_", regardless of what the
// integration error says.  This is used by advance_fixed
template <class System>
void SolverInternal<System>::step_to(System& system, double time_max_) {
  set_time_max(time_max_);
  replay_step(system); // restore this step's recorded state on a replay pass
  setup_dydt_in(system);
  stepper_step(system, time, time_max - time, y, yerr, dydt_in, dydt_out);
  save_dydt_out_as_in();
  record_ode_step(system);

  time = time_max;
  prev_times.push_back(time);
}

template <class System>
void SolverInternal<System>::resize(size_t size_) {
  y.resize(size_);
  yerr.resize(size_);
  dydt_in.resize(size_);
  dydt_out.resize(size_);
  stepper.resize(size_);
  rodas_stepper.resize(size_);
}

template <class System>
void SolverInternal<System>::setup_dydt_in(System& system) {
  if (stepper_can_use_dydt_in() && !dydt_in_is_clean) {
    // TODO: Not clear that this is the right thing here; should just
    // be able to look up the correct dydt rates because we've already
    // set state?
    ode::derivs(system, y, dydt_in, time);
    dydt_in_is_clean = true;
  }
}

template <class System>
void SolverInternal<System>::save_dydt_out_as_in() {
  if (stepper_first_same_as_last()) {
    dydt_in = dydt_out;
    dydt_in_is_clean = true;
  } else {
    dydt_in_is_clean = false;
  }
}

template <typename System>
void SolverInternal<System>::set_time(double t) {
  const int ulp = 2; // units in the last place (accuracy)
  if (prev_times.size() > 0 &&
      !util::almost_equal(prev_times.back(), t, ulp))
  {
    util::stop("Time does not match previous (delta = " +
               util::to_string(prev_times.back() - t) +
               "). Reset solver first.");
  }
  time = t;
  if (prev_times.empty()) { // only if first time (avoids duplicate times)
    prev_times.push_back(time);
  }
}

template <class System>
void SolverInternal<System>::set_time_max(double time_max_) {
  if (!util::is_finite(time_max_)) {
    util::stop("time_max must be finite!");
  }
  if (time_max_ < time) {
    util::stop("time_max must be greater than (or equal to) current time");
  }
  time_max = time_max_;
}

}
}

#endif