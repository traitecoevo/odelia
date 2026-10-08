// -*-c++-*-
#ifndef ODELIA_ODE_SOLVER_INTERNAL_HPP_
#define ODELIA_ODE_SOLVER_INTERNAL_HPP_

#include <array>
#include <odelia/ode_interface.hpp>
#include <odelia/ode_control.hpp>
#include <odelia/ode_step_rkck.hpp>
#include <odelia/ode_step_rodas.hpp>
#include <odelia/ode_step_dopri.hpp>

#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <span>
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
  std::vector<double> get_times() const;
  std::vector<double> get_step_sizes() const;

  void advance_adaptive(System &system, double time_max_);
  void advance_fixed(System& system, const std::vector<double>& times);
  void advance_euler(System& system, const std::vector<double>& times);

  void step(System& system);

  // The adjoint of one step, from the state that step started at, for several
  // seeds at once. RKCK only: neither the Rosenbrock nor the Dormand-Prince
  // stepper records its stages, so neither carries a reverse
  // counterpart.
  void step_adjoint(active_system<System>& active,
                    const typename Step<System>::solved_row& solved, double time,
                    double step_size, const state_type& y,
                    const adjoint_rows& lambda_out, adjoint_rows& lambda_in,
                    adjoint_rows& parameter_adjoint) {
    if (method == Method::rodas) {
      util::stop("method='rodas' has no adjoint; use method='rkck'.");
    }
    if (method == Method::dopri) {
      util::stop("method='dopri' has no adjoint; use method='rkck'.");
    }
    // The System can be a different width from the one the forward pass left,
    // because a caller sweeping a range narrows it between ranges. Every seed
    // carries that same width, and the stage buffers are sized to it here. A sweep
    // is the end of the solver's forward state either way.
    resize(lambda_out.width());
    stepper.step_adjoint(active, solved, time, step_size, y,
                         lambda_out, lambda_in, parameter_adjoint);
  }

  // Rate evaluations recorded since the count was last cleared.
  std::size_t recorded_rates() const { return stepper.recorded_rates; }
  void clear_recorded_rates() { stepper.recorded_rates = 0; }


  // Keep the state at each accepted step as well as the time and the size. The
  // caller's decision: a run whose gradient will be taken needs the states, and a
  // run that is only integrating does not.
  //
  // Set before the first step. The row the run starts from is already in the
  // schedule by then (reset() put it there), so turning this on fills that row's
  // state from the state the solver holds, and the record covers every step.
  void set_keep_states(bool keep) {
    keep_states_ = keep;
    if (keep && prev_schedule.size() == 1 && prev_steps.empty()) {
      step_record<System> row{prev_schedule.front(), state_type()};
      row.state.assign(y.begin(), y.end());
      prev_steps.push_back(std::move(row));
    }
  }
  // The record itself, which is what a sweep reads. One row per accepted step,
  // carrying the time, the size that reached it and the state there -- so a
  // caller cannot pair one run's state with another run's size, and cannot be
  // handed a time without the state it belongs to.
  //
  // Refused where the run was not asked to keep states, because a row without
  // one is not a row a sweep can use and returning it would move the check to
  // every caller.
  std::span<const step_record<System>> recording() const {
    if (!keep_states_) {
      util::stop("recording(): this run was not asked to keep its states, so "
                 "there is no record to sweep -- set_keep_states(true) before "
                 "the run");
    }
    if (prev_steps.size() != prev_schedule.size()) {
      util::stop("recording(): keep_states was turned on after the run had "
                 "stepped, so the record does not cover every step -- "
                 "set_keep_states(true) before the first step");
    }
    return {prev_steps.data(), prev_steps.size()};
  }

  // The steps this run took: each time it reached and the size that reached it,
  // paired as the run paired them. Available whether or not states were kept,
  // because what the run decided is not what it held.
  //
  // Insertion rows are left out, and this is the one place that decides so. An
  // insertion shares its time with the row below it, so a caller pinning these
  // times would be handed one twice -- and the two callers that pin them
  // (a schedule written back into Parameters, and a replay driven interval by
  // interval) both apply their own insertions. get_times and get_step_sizes are
  // this list split in two rather than two more walks that could disagree.
  std::vector<instruction> schedule() const {
    std::vector<instruction> ret;
    ret.reserve(prev_schedule.size());
    for (const instruction& s : prev_schedule) {
      if (!s.insertion) {
        ret.push_back(s);
      }
    }
    return ret;
  }

  // One accepted step, into the record: the time it reached, the size that
  // reached it, and the state there where the run was asked to keep states.
  // `subdivided` marks a pinned interval crossed in more than one sub-step.
  void push_step(System& system, double time_, double step_size,
                 bool subdivided = false);
  // The insertion the caller just applied, as a row of its own: it holds the state
  // the map produced, at the time the row below it holds. Only schedule() has to
  // know that two rows share a time, and it drops these.
  void push_insertion(System& system);
  void step(System& system, double time_max_);
  void step_to(System& system, double time_max_);
  // `reached` is the time a recording says this step ended at; NaN accumulates.
  // `replay`, where a recording is being walked, is that step's own solved row.
  void step_by(System& system, double step_size, double reached,
               const typename Step<System>::solved_row* replay = nullptr);
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

  // Dense output: the state at any t inside the last accepted adaptive
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
  // The step index the stages are addressed by is this object's own count of
  // accepted steps, which is the step about to be taken -- read here rather than
  // passed, so no caller can disagree with it.
  //
  // `replay` is the row an earlier run wrote, where a caller is replaying one:
  // the stages then load it instead of solving again, and nothing is cleared
  // because nothing is written.
  void stepper_step(System& system, double time_, double step_size,
                    state_type& y_, state_type& yerr_,
                    const state_type& dydt_in_, state_type& dydt_out_,
                    const typename Step<System>::solved_row* replay = nullptr) {
    if (method == Method::dopri) {
      if (replay != nullptr) {
        util::stop("method='dopri' cannot replay what a run solved for: the "
                   "Dormand-Prince stepper keeps no per-stage row; use "
                   "method='rkck'.");
      }
      dopri_stepper.step(system, time_, step_size, y_, yerr_, dydt_in_, dydt_out_);
    } else if (method == Method::rodas) {
      if (replay != nullptr) {
        // RODAS takes no row at all, so it has nowhere to put one. Said here
        // rather than dropped silently, which would replay a run's choices by
        // re-deriving them.
        util::stop("method='rodas' cannot replay what a run solved for: the "
                   "Rosenbrock stepper keeps no per-stage row; use method='rkck'.");
      }
      if constexpr (RodasStep<System>::supported) {
        rodas_stepper.step(system, time_, step_size, y_, yerr_, dydt_in_,
                           dydt_out_);
      } else {
        // RODAS is unavailable for this System: it has neither an
        // ode_jacobian() hook nor a rebind_from() hook for the AD Jacobian, or its
        // scalar type is itself active.
        util::stop("method='rodas' is not available for this system/scalar type "
                   "(needs an ode_jacobian() hook, or a rebind_from() hook with a "
                   "non-active scalar); use method='rkck'.");
      }
    } else {
      // Into scratch, because a step that is rejected and retried writes here twice
      // and only the one that is accepted is committed -- which is what makes
      // "a rejected attempt writes the same slot as its retry" nothing anyone has
      // to arrange.
      if (replay != nullptr) {
        stepper.step(system, *replay, time_, step_size, y_, yerr_, dydt_in_,
                     dydt_out_);
      } else {
        for (solved_values_t<System>& row : solved_scratch_) { row = {}; }
        stepper.step(system, solved_scratch_, time_, step_size, y_, yerr_,
                     dydt_in_, dydt_out_);
      }
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
  // Each accepted step: the time it reached and the size it took, both recorded
  // (see instruction). One entry per accepted step, the first being the state the run started from,
  // which no step reached and which therefore has no size.
  //
  // The state is kept BESIDE the size that reached it because the two are one
  // record: in separate stores a walk could pair a state with a size from a
  // different run, and nothing would say so.
  std::vector<step_record<System>> prev_steps;
  // What every run keeps, recording or not: the time each step reached, the
  // size that reached it, and where the state widened. Separate from the
  // records so a run that is only integrating stores 24 bytes a step rather
  // than a record with a state and a row of solved values. The two are written
  // together, row for row, and recording() refuses if they differ.
  std::vector<instruction> prev_schedule;
  // Whether to keep the states. The caller's: a run whose gradient will be taken
  // needs them and a run that is only integrating does not. The same flag decides
  // whether what a step solves for is kept, because those are one recording.
  bool keep_states_ = false;

  // The step being attempted writes what it solves for here, and an accepted step
  // moves it onto its row.
  typename Step<System>::solved_row solved_scratch_;

  state_type y;        // Vector of current system state
  state_type yerr;     // Vector of error estimates
  state_type dydt_in;  // Vector of dydt at beginning of step
  state_type dydt_out; // Vector of dydt during step

  bool dydt_in_is_clean;
};

template <class System>
SolverInternal<System>::SolverInternal(System &system, OdeControl control_,
                                       Method method_)
  : control(control_), method(method_) {
  reset(system);
}

// NOTE: This resets *everything* to basically a recreated object.
template <class System>
void SolverInternal<System>::reset(System& system) {
  prev_steps.clear();
  prev_schedule.clear();
  step_size_last = control.step_size_initial;
  time_max = std::numeric_limits<double>::infinity();
  in_step = false;
  have_prev = false;
  set_state_from_system(system);
}

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
  // The state the run starts from, into the record set_time just opened. Kept
  // here rather than by the caller because this is the one place that knows the
  // System has been read.
  if (keep_states_ && prev_steps.size() == 1 && prev_steps.back().state.empty()) {
    prev_steps.back().state.assign(y.begin(), y.end());
  }
}

// One accepted step, recorded. The state comes off the System rather than out of
// `y`, so a System that reaches a state by a route of its own is recorded at the
// state it holds.
template <class System>
void SolverInternal<System>::push_step(System& system, double time_,
                                       double step_size, bool subdivided) {
  prev_schedule.push_back({time_, step_size});
  if (keep_states_) {
    step_record<System> record{{time_, step_size}, state_type()};
    record.state.resize(system.ode_size());
    system.ode_state(record.state.begin());
    record.solved = std::move(solved_scratch_);
    record.subdivided = subdivided;
    prev_steps.push_back(std::move(record));
  }
}

// The row goes in on every run, because where a run widened is what it decided;
// its state is filled only where states are kept, because nothing but a sweep
// reads one.
template <class System>
void SolverInternal<System>::push_insertion(System& system) {
  if (prev_schedule.empty()) {
    util::stop("push_insertion: no recorded step for an insertion to follow");
  }
  const instruction row{prev_schedule.back().time,
                        std::numeric_limits<double>::quiet_NaN(), true};
  prev_schedule.push_back(row);
  if (keep_states_) {
    step_record<System> record{row, state_type()};
    record.state.resize(system.ode_size());
    system.ode_state(record.state.begin());
    prev_steps.push_back(std::move(record));
  }
}

template <class System>
std::vector<double> SolverInternal<System>::get_times() const {
  const std::vector<instruction> steps = schedule();
  std::vector<double> ret;
  ret.reserve(steps.size());
  for (const instruction& s : steps) {
    ret.push_back(s.time);
  }
  return ret;
}

// The size of the step that reached each recorded time; NaN for the initial
// time, which no step reached.
template <class System>
std::vector<double> SolverInternal<System>::get_step_sizes() const {
  const std::vector<instruction> steps = schedule();
  std::vector<double> ret;
  ret.reserve(steps.size());
  for (const instruction& s : steps) {
    ret.push_back(s.step_size);
  }
  return ret;
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
// size adjustments will be done.
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
  push_step(system, time, h);
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
    // usually lands inside the domain.
    //
    //   1. A stage throws util::DomainError, which is how a model reports an
    //      out-of-domain state.
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
          // with the system and y disagreeing.
          //
          // Deliberately not done on an accuracy rejection: there the system sits
          // on the completed step's final state and the retry's stage 2 overwrites
          // it before anything reads it, so a state-set -- which can be costly for
          // a System -- would buy nothing.
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
      time_prev = time_orig;
      have_prev = true;
      save_dydt_out_as_in();
      push_step(system, time, step_size);
      in_step = false;
      return; // This exits the infinite loop.
    }
  }
}

// One adaptive step that will not pass time_max_: the single-step form of
// advance_adaptive(), for a caller that drives the integration itself and may
// change the state between steps. An infinite time_max_ removes the bound,
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
// ways step() refuses one: a stage throws util::DomainError, or the completed
// step lands on a state ode_state_valid() refuses. step() answers either by
// rejecting the step and retrying it smaller. Here the endpoint is given by the
// caller and cannot be moved, so the interval is subdivided instead -- shrink the
// sub-step and walk to the same endpoint in several. The endpoint is still hit
// exactly, so the times the caller records are unchanged, and an interval that
// raises neither objection takes exactly one step. Without it, a pinned run of a
// System that refuses states routinely could not complete.
//
// A subdivided interval is still one row of the record, holding the last
// sub-step's solved values, and the row says so (step_record::subdivided): a
// sweep or a replay refuses it rather than treating the interval as one step.
template <class System>
void SolverInternal<System>::step_to(System& system, double time_max_) {
  set_time_max(time_max_);
  in_step = true;
  have_prev = false; // a pinned step may be subdivided; no single interpolant
  // The interval, which is what the recording holds however many sub-steps the
  // retry below takes to cross it: a replay reproduces the caller's times from
  // this row, and a shrunken sub-step is a detail of how this run got there.
  const double step_size = time_max - time;
  setup_dydt_in(system);

  // The sub-step, named apart from the interval above because the retry
  // shrinks it. Starts as the whole interval, so the common case is one step.
  double sub_step_size = time_max - time;
  // Held across iterations so a retry reuses the buffer rather than allocating.
  state_type y_orig;
  // Sub-steps accepted; more than one marks the row (step_record::subdivided).
  size_t sub_steps = 0;

  while (true) {
    // Take the endpoint from time_max rather than accumulating sub_step_size,
    // so the caller's time is reproduced bit-for-bit however the interval was
    // cut. A zero-length interval lands here immediately and steps once.
    const bool final_sub_step = !(time + sub_step_size < time_max);
    const double time_next = final_sub_step ? time_max : time + sub_step_size;

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
    // Both ways of refusing a state apply here, as on the adaptive path;
    // honouring only the throw would leave the predicate silently unenforced
    // whenever the integration happens to be pinned.
    if (!invalid && !state_valid(system, y)) {
      invalid = true;
      invalid_reason = "ode_state_valid() refused the state after the step";
    }

    if (!invalid) {
      save_dydt_out_as_in();
      time = time_next;
      ++sub_steps;
      if (final_sub_step) {
        break;
      }
      continue;
    }

    // Undo the failed sub-step. The system is left holding whichever stage threw,
    // so put it back on the restored state explicitly -- otherwise a give-up below
    // would exit with the system and y disagreeing.
    y = y_orig;
    internal::set_ode_state(system, y, time);

    // Shrink the same way an invalid step shrinks on the adaptive path, and stop
    // at the same floor: one rule for how hard the solver retries a state a system
    // refused, and one knob (step_size_min) controlling it.
    const double step_size_next = control.reject_step(sub_step_size);
    if (!(step_size_next < sub_step_size) || !(time + step_size_next > time)) {
      // The smallest sub-step we are allowed to take still leaves the domain (or
      // is too small to advance the clock at all). Name the reason, because it came
      // from the system and is the only description of what is actually wrong.
      util::stop("Cannot leave an invalid state at t = " +
                 util::format_double(time) + " stepping to " +
                 util::format_double(time_max) + ": " + invalid_reason +
                 " (sub-step size " + util::format_double(sub_step_size) +
                 " is already at the minimum)");
    }
    sub_step_size = step_size_next;
    ++n_rejections;
  }

  time = time_max;
  push_step(system, time, step_size, sub_steps > 1);
  in_step = false;
}

// This takes a step of the given size, regardless of what the integration error
// says. advance_recorded() drives a replay through it, one recorded step each.
template <class System>
void SolverInternal<System>::step_by(System& system, double step_size,
                                     double reached,
                                     const typename Step<System>::solved_row*
                                       replay) {
  if (!util::is_finite(step_size)) {
    util::stop("step_size must be finite!");
  }
  if (step_size < 0.0) {
    util::stop("step_size must be greater than (or equal to) zero");
  }
  in_step = true;
  have_prev = false; // a fixed step leaves no adaptive step to interpolate in
  setup_dydt_in(system);
  stepper_step(system, time, step_size, y, yerr, dydt_in, dydt_out, replay);
  save_dydt_out_as_in();

  // The time the run reached, where a recording says what it was, rather than
  // this time plus this size (see instruction).
  time = util::is_finite(reached) ? reached : time + step_size;
  time_max = time;
  push_step(system, time, step_size);
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
  // Only the stepper that will run. `method` is fixed at construction and there
  // is no setter, so the others' scratch is never read -- and the Rosenbrock
  // scratch is two size x size matrices, which a sweep would otherwise zero once
  // per recorded step.
  if (method == Method::rodas) {
    rodas_stepper.resize(size_);
  } else if (method == Method::dopri) {
    dopri_stepper.resize(size_);
  } else {
    stepper.resize(size_);
  }
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
  if (prev_schedule.size() > 0 &&
      !util::almost_equal(prev_schedule.back().time, t, ulp))
  {
    util::stop("Time does not match previous (delta = " +
               util::format_double(prev_schedule.back().time - t) +
               "). Reset solver first.");
  }
  time = t;
  if (prev_schedule.empty()) { // only if first time (avoids duplicate times)
    // No step reached the initial time, so it records no size. The state is
    // recorded by set_state_from_system, which calls this and then holds it.
    const instruction row{time, std::numeric_limits<double>::quiet_NaN()};
    prev_schedule.push_back(row);
    if (keep_states_) {
      prev_steps.push_back(step_record<System>{row, state_type()});
    }
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