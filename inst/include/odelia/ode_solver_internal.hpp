// -*-c++-*-
#ifndef ODELIA_ODE_SOLVER_INTERNAL_HPP_
#define ODELIA_ODE_SOLVER_INTERNAL_HPP_

#include <odelia/ode_interface.hpp>
#include <odelia/ode_control.hpp>
#include <odelia/ode_step_rkck.hpp>
#include <odelia/ode_step_rodas.hpp>
#include <odelia/ode_step_dopri.hpp>

#include <limits>
#include <string>
#include <vector>
#include <cstddef>

namespace odelia {
namespace ode {

// Integration method: the explicit Cash-Karp RKCK 4(5) stepper (default), the
// implicit RODAS4(3) Rosenbrock stepper for stiff systems, or the explicit
// Dormand-Prince 5(4) stepper with order-4 dense output.
enum class Method { rkck, rodas, dopri };

template <class System>
class SolverInternal {
public:
  // Extract scalar type from System using traits
  using value_type = typename System::value_type;
  using state_type = std::vector<value_type>;

  // The system is taken by mutable reference throughout: reading its rates may
  // require it to compute them for the state it currently holds. See
  // set_state_from_system.
  SolverInternal(System &system, OdeControl control_,
                 Method method_ = Method::rkck);
  void reset(System& system);
  void set_state_from_system(System& system);

  state_type get_state() const {return y;}
  double get_time() const {return time;}
  std::vector<double> get_times() const {return prev_times;}

  void advance_adaptive(System &system, double time_max_);
  void advance_fixed(System& system, const std::vector<double>& times);
  void advance_euler(System& system, const std::vector<double>& times);

  void step(System& system);
  void step(System& system, double time_max_);
  void step_to(System& system, double time_max_);
  void step_euler(System& system, double time_max_);

  void set_time_max(double time_max_);

  // The step the controller will try next: the size of the last accepted step
  // as adjusted by its error estimate, or the control's initial size after a
  // reset. A consumer that re-seeds the state through set_state(), which resets
  // this, may put a step it knows to be good back with set_step_size().
  double get_step_size() const { return step_size_last; }
  void set_step_size(double h) {
    if (!util::is_finite(h) || h <= 0.0) {
      util::stop("step size must be positive and finite");
    }
    step_size_last = h;
  }
  // Attempts the adaptive and pinned paths rejected and retried smaller, since
  // construction. A diagnostic: with the system's own count of right-hand-side
  // evaluations it says what a solve cost and why.
  size_t get_n_rejections() const { return n_rejections; }
  // True while an attempt is in progress, and therefore afterwards if one was
  // abandoned by an exception that escaped step(): y and the system then hold
  // a half-finished attempt and must be re-seeded (reset()) before stepping
  // again. A consumer whose right-hand side can raise -- an R callback -- reads
  // this rather than guessing from where the error came from.
  bool mid_step() const { return in_step; }

  // Dense output (#24): the state at any t inside the last accepted adaptive
  // step, at no further evaluation. Under Dormand-Prince it is the method's
  // own order-4 continuous extension, from the step's stages, with an error of
  // the step's own order. Under the other steppers it is cubic Hermite on the
  // step's endpoints and the derivatives there: exact for a cubic, otherwise
  // one order below the step (measured on Lorenz at 40-100x the integration
  // error), so a consumer that wants dense output within tolerance uses
  // Method::dopri. Refuses a t outside [previous time, current time], and is
  // unavailable until a step has been accepted since the last reset.
  bool can_interpolate() const { return have_prev; }
  double previous_time() const { return time_prev; }
  void interpolate(double t, state_type& out) const;

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
    if (method == Method::dopri) {
      dopri_stepper.step(system, time_, step_size, y_, yerr_, dydt_in_, dydt_out_);
    } else if (method == Method::rodas) {
      if constexpr (RodasStep<System>::supported) {
        rodas_stepper.step(system, time_, step_size, y_, yerr_, dydt_in_,
                           dydt_out_);
      } else {
        // RODAS is unavailable for this System: it has neither an
        // ode_jacobian() hook nor a rebind() hook for the AD Jacobian, or its
        // scalar type is itself active (nested tangent-over-adjoint is not yet
        // wired up -- see issue #36).
        util::stop("method='rodas' is not available for this system/scalar type "
                   "(needs an ode_jacobian() hook, or a rebind() hook with a "
                   "non-active scalar); use method='rkck'.");
      }
    } else {
      stepper.step(system, time_, step_size, y_, yerr_, dydt_in_, dydt_out_);
    }
  }
  size_t stepper_order() const {
    switch (method) {
    case Method::rodas: return rodas_stepper.order();
    case Method::dopri: return dopri_stepper.order();
    default: return stepper.order();
    }
  }
  bool stepper_can_use_dydt_in() const {
    switch (method) {
    case Method::rodas: return RodasStep<System>::can_use_dydt_in;
    case Method::dopri: return DopriStep<System>::can_use_dydt_in;
    default: return Step<System>::can_use_dydt_in;
    }
  }
  bool stepper_first_same_as_last() const {
    switch (method) {
    case Method::rodas: return RodasStep<System>::first_same_as_last;
    case Method::dopri: return DopriStep<System>::first_same_as_last;
    default: return Step<System>::first_same_as_last;
    }
  }

  OdeControl control;
  Method method;
  Step<System> stepper;
  RodasStep<System> rodas_stepper;
  DopriStep<System> dopri_stepper;

  double step_size_last; // Size of last successful step (or suggestion)
  size_t n_rejections = 0; // Rejected attempts, cumulative (not reset)
  bool in_step = false;    // An attempt is in progress (see mid_step())

  // The last accepted adaptive step: its start state and derivative, for the
  // interpolant and for restoring y on a rejected attempt.
  state_type y_prev;
  state_type dydt_prev;
  double time_prev = 0.0;
  bool have_prev = false;

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
SolverInternal<System>::SolverInternal(System &system, OdeControl control_,
                                       Method method_)
  : control(control_), method(method_) {
  reset(system);
}

// NOTE: This resets *everything* to basically a recreated object.
template <class System>
void SolverInternal<System>::reset(System& system) {
  prev_times.clear();
  step_size_last = control.step_size_initial;
  time_max = std::numeric_limits<double>::infinity();
  in_step = false;
  have_prev = false;
  set_state_from_system(system);
}

// saving ode steps during adaptive solve
template <typename System>
typename std::enable_if<has_cache<System>::value, void>::type
cache(System& system) {
  system.cache_ode_step();
}

template <typename System>
typename std::enable_if<!has_cache<System>::value, void>::type
cache(System& system) {}

// During mutant run, load ode history
// the history is a vector of 6 states for env, needed to make a 
// full RK step. Called as part of `step_to`
template <typename System>
typename std::enable_if<has_cache<System>::value, void>::type
load(System& system) {
  system.load_ode_step();
}

// During resident run, no cache loaded, proceed as normal
template <typename System>
typename std::enable_if<!has_cache<System>::value, void>::type
load(System& system) {}

// Seed y and dydt_in from whatever state the system currently holds. The system
// is mutable because `ode_rates` is allowed to compute: a system that reaches a
// state by a route of its own (widening it, reloading it) can then hand back the
// derivative *of that state* rather than a cached one belonging to an earlier
// one. Marking dydt_in clean here is only sound because of that -- with a const
// system the rates were whatever the system last happened to store, and under
// first-same-as-last they became k1 of the next step.
template <class System>
void SolverInternal<System>::set_state_from_system(System& system) {
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
  // changes 'y'). Kept as a member: it is also the start of the step the
  // interpolant reads, once the step is accepted.
  y_prev = y;
  const state_type& y_orig = y_prev;
  const size_t size = y.size();

  in_step = true;
  // Compute the derivatives at the beginning.
  setup_dydt_in(system);
  dydt_prev = dydt_in;

  while (true) {
    // Does this appear to be the last step before reaching `time_max`?
    const bool final_step = step_size > time_remaining;
    if (final_step) {
      step_size = time_remaining;
    }

    // Beyond being inaccurate, a step can be *invalid* in two ways, and both are
    // rejections rather than failures: y_orig is right here, and a smaller step
    // usually lands inside the domain (#55).
    //
    //   1. A stage throws util::DomainError. This is how a model normally reports
    //      an out-of-domain state, and until now such a throw escaped this
    //      function and killed the whole solve.
    //   2. The completed step lands on a state the system's optional
    //      ode_state_valid() refuses.
    //
    // Only DomainError is caught. Anything else -- util::stop(), std::bad_alloc,
    // a logic error -- still propagates, so a bug stays a bug instead of becoming
    // "Cannot achieve the desired accuracy".
    bool invalid = false;
    std::string invalid_reason;
    try {
      stepper_step(system, time, step_size, y, yerr, dydt_in, dydt_out);
    } catch (const util::DomainError& e) {
      invalid = true;
      invalid_reason = e.what();
    }

    double step_size_next;
    if (invalid) {
      // yerr and dydt_out were never completed, so there is no error estimate to
      // form: reject on the strength of the throw alone.
      step_size_next = control.reject_step(step_size);
    } else {
      step_size_next =
        control.adjust_step_size(size, stepper_order(), step_size,
			         y, yerr, dydt_out);
      if (!state_valid(system, y)) {
        invalid = true;
        invalid_reason = "ode_state_valid() refused the state after the step";
        // Overrides whatever the error estimate concluded, including "accept".
        step_size_next = control.reject_step(step_size);
      }
    }

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
        ++n_rejections;
        if (invalid) {
          // Put the system back on the restored state explicitly. After a caught
          // DomainError it is left holding whichever intermediate stage threw, and
          // if the retry goes on to raise at the minimum step size we would exit
          // with the system and y disagreeing -- the pattern behind the stale-state
          // bugs (plant#585, plant#589).
          //
          // Deliberately not done on an accuracy rejection: there the system sits
          // on the completed step's final state, the retry's stage 2 overwrites it
          // before anything reads it, and that has always been the behaviour. Doing
          // it unconditionally would add a state-set -- for plant, an environment
          // rebuild -- to every rejected step, for systems that gain nothing from
          // this feature.
          internal::set_ode_state(system, y, time);
        }
      } else {
      	// We've reached limits of machine accuracy in differences of
      	// step sizes or time (or both).
        if (invalid) {
          // Not an accuracy problem: the smallest step we are allowed to take
          // still leaves the domain. Name the reason, because it came from the
          // system and is the only description of what is actually wrong.
          util::stop("Cannot leave an invalid state at t = " +
                     util::format_double(time_orig) + ": " + invalid_reason +
                     " (step size " + util::format_double(step_size) +
                     " is already at the minimum)");
        }
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
      time_prev = time_orig;
      have_prev = true;
      save_dydt_out_as_in();
      cache(system);
      in_step = false;
      return; // This exits the infinite loop.
    }
  }
}

// One adaptive step that will not pass time_max_: the single-step form of
// advance_adaptive(), for a caller that drives the integration itself and may
// change the state between steps (#62). An infinite time_max_ removes the bound,
// which is what a solver has after reset(). Stepping from time_max_ itself is
// refused: a zero-length step is not a step, and the implicit stepper divides by
// h.
template <class System>
void SolverInternal<System>::step(System& system, double time_max_) {
  if (util::is_finite(time_max_)) {
    set_time_max(time_max_);
    if (!(time < time_max)) {
      util::stop("step(): already at time_max = " + util::format_double(time_max));
    }
  } else {
    time_max = std::numeric_limits<double>::infinity();
  }
  step(system);
}

// This takes a step up to time "time_max_", regardless of what the
// integration error says.  This is used by advance_fixed
//
// The step is not error-controlled, but it can still be *invalid*, in the two
// ways #55 defines: a stage throws util::DomainError, or the completed step lands
// on a state ode_state_valid() refuses. step() answers either by rejecting the
// step and retrying it smaller. Here the endpoint is given by the caller and
// cannot be moved, so the interval is subdivided instead -- shrink the sub-step
// and walk to the same endpoint in several. The endpoint is still hit exactly, so
// the times the caller records are unchanged, and an interval that raises neither
// objection takes exactly one step, as before.
//
// Without this, the rejection added in #55 was reachable only from the adaptive
// path: advance_fixed called the stepper bare, so the first throw killed the
// solve. That made a whole class of run impossible rather than slow -- plant's
// mutant replay pins the stepper to a resident's recorded times, and its TF24
// model throws here as a matter of routine (~480 rejections in a resident run
// that goes on to complete), so a replay was near-certain to meet one
// (plant#642).
//
// One caveat for systems that cache per-stage data (`cache(system, rk_step)`):
// the stage indices restart at 0 on each sub-step, so a subdivided interval
// leaves the system holding the *last* sub-step's stages rather than stages
// spanning the whole interval. Consumers that record such a cache for later
// replay get a coarser record of a subdivided step than of a plain one.
template <class System>
void SolverInternal<System>::step_to(System& system, double time_max_) {
  set_time_max(time_max_);
  in_step = true;
  have_prev = false; // a pinned step may be subdivided; no single interpolant
  load(system); // option to load pre-calculated states in mutant runs
  setup_dydt_in(system);

  // Sub-step size. Starts as the whole interval, so the common case is one step.
  double step_size = time_max - time;
  // Held across iterations so a retry reuses the buffer rather than allocating.
  state_type y_orig;

  while (true) {
    // Take the endpoint from time_max rather than accumulating step_size, so the
    // caller's time is reproduced bit-for-bit however the interval was cut. A
    // zero-length interval lands here immediately and steps once, as before.
    const bool final_sub_step = !(time + step_size < time_max);
    const double time_next = final_sub_step ? time_max : time + step_size;

    y_orig = y;
    bool invalid = false;
    std::string invalid_reason;
    try {
      // dydt_in is read, not written, so only y needs saving to retry.
      stepper_step(system, time, time_next - time, y, yerr, dydt_in, dydt_out);
    } catch (const util::DomainError& e) {
      invalid = true;
      invalid_reason = e.what();
    }
    // The other half of the #55 contract. Both ways of refusing a state apply
    // here for the same reason they apply on the adaptive path; honouring only
    // the throw would leave the predicate silently unenforced whenever the
    // integration happens to be pinned.
    if (!invalid && !state_valid(system, y)) {
      invalid = true;
      invalid_reason = "ode_state_valid() refused the state after the step";
    }

    if (!invalid) {
      save_dydt_out_as_in();
      time = time_next;
      if (final_sub_step) {
        break;
      }
      continue;
    }

    // Undo the failed sub-step. The system is left holding whichever stage threw,
    // so put it back on the restored state explicitly -- otherwise a give-up below
    // would exit with the system and y disagreeing, the pattern behind the
    // stale-state bugs (plant#585, plant#589).
    y = y_orig;
    internal::set_ode_state(system, y, time);

    // Shrink the same way an invalid step shrinks on the adaptive path, and stop
    // at the same floor: one rule for how hard the solver retries a state a system
    // refused, and one knob (step_size_min) controlling it.
    const double step_size_next = control.reject_step(step_size);
    if (!(step_size_next < step_size) || !(time + step_size_next > time)) {
      // The smallest sub-step we are allowed to take still leaves the domain (or
      // is too small to advance the clock at all). Name the reason, because it came
      // from the system and is the only description of what is actually wrong.
      util::stop("Cannot leave an invalid state at t = " +
                 util::format_double(time) + " stepping to " +
                 util::format_double(time_max) + ": " + invalid_reason +
                 " (sub-step size " + util::format_double(step_size) +
                 " is already at the minimum)");
    }
    step_size = step_size_next;
    ++n_rejections;
  }

  cache(system);

  time = time_max;
  prev_times.push_back(time);
  in_step = false;
}

template <class System>
void SolverInternal<System>::interpolate(double t, state_type& out) const {
  if (!have_prev) {
    util::stop("interpolate(): no accepted step to interpolate within");
  }
  if (!(t >= time_prev && t <= time)) {
    util::stop("interpolate(): t = " + util::format_double(t) +
               " is outside the last step [" + util::format_double(time_prev) +
               ", " + util::format_double(time) + "]");
  }
  const size_t n = y.size();
  out.resize(n);
  const double h = time - time_prev;
  if (!(h > 0.0)) {
    out = y;
    return;
  }
  const double s = (t - time_prev) / h;
  if (method == Method::dopri) {
    dopri_stepper.dense(s, y, out);
    return;
  }
  const double s2 = s * s, s3 = s2 * s;
  const double h00 = 2.0 * s3 - 3.0 * s2 + 1.0;
  const double h10 = (s3 - 2.0 * s2 + s) * h;
  const double h01 = -2.0 * s3 + 3.0 * s2;
  const double h11 = (s3 - s2) * h;
  // dydt_in is f at the current y: FSAL carried it, or the next step's
  // setup will compute it. Only the carried case can be read here.
  const state_type& dydt_now = dydt_in_is_clean ? dydt_in : dydt_out;
  for (size_t i = 0; i < n; ++i) {
    out[i] = h00 * y_prev[i] + h10 * dydt_prev[i] + h01 * y[i] + h11 * dydt_now[i];
  }
}

template <class System>
void SolverInternal<System>::resize(size_t size_) {
  y.resize(size_);
  yerr.resize(size_);
  dydt_in.resize(size_);
  dydt_out.resize(size_);
  stepper.resize(size_);
  rodas_stepper.resize(size_);
  dopri_stepper.resize(size_);
}

template <class System>
void SolverInternal<System>::setup_dydt_in(System& system) {
  if (stepper_can_use_dydt_in() && !dydt_in_is_clean) {
    // The full derivs() -- set the state, then read the rates -- rather than
    // reading the rates alone. y is the solver's own state and need not be the
    // one the system currently holds, so the state has to be re-established
    // before the rates mean anything.
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