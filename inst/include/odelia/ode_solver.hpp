#ifndef ODELIA_ODE_SOLVER_HPP_
#define ODELIA_ODE_SOLVER_HPP_

#include <odelia/ode_solver_internal.hpp>
#include <odelia/adjoint.hpp>
#include <algorithm>
#include <utility>
#include <XAD/XAD.hpp>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <optional>
#include <span>
#include <vector>

namespace odelia {
namespace ode {

// Solver<System> owns a System and a SolverInternal and drives them: forward
// over a set of times (advance_adaptive, advance_fixed, advance_schedule, replay),
// and backward over what a forward run recorded (sweep_adjoint). The System is
// a copy of the one the Solver was built from, which is left as it was; the
// copy is the one that steps, and get_system_ref() is how a caller reaches it.
//
// A gradient of a run, in full:
//
//   Solver<Lorenz> s(system, control);
//   s.set_keep_states(true);              // keep one step_record per step
//   s.advance_adaptive({0.0, 10.0});      // the ordinary solve
//   adjoint_batch lambda = adjoint_batch::one({1.0, 0.0, 0.0});  // d(x_final)
//   adjoint_batch dp(1, n_parameters);     // one per output, zeroed
//   s.sweep_adjoint(lambda, dp);
//   // dp[0][j] = d x_final / d parameter_j; lambda[0][i] = d x_final / d y_i(0)
//
// Where to read, in order: step_record and the Sweepable concept
// (ode_interface.hpp) for what a run keeps and what a System must have;
// push_step and push_insertion (ode_solver_internal.hpp) for the two places a
// row is written; state_and_parameter_adjoints (adjoint.hpp) for the transpose
// of one map; Step::step_adjoint (ode_step_rkck.hpp) for the transpose of one
// step, which DopriStep::step_adjoint (ode_step_dopri.hpp) repeats on its own
// tableau; sweep_adjoint below for the loop over rows and insertions. Under
// method = rkck or dopri a run records, replays and sweeps; rodas refuses a sweep
// or a replay.
//
// The same recording is differentiated forward, for few inputs and many
// outputs and under any of the three steppers, by sweep_tangent: see the free
// sweep_tangent below, and the article "Forward mode" for a worked example.
//
// `collect` (default on) keeps a copy of the whole System at every step in
// `history`, for callers reading a trajectory from R. A sweep does not read it
// and does not need it; set_collect(false) where the copies cost.

template <std::size_t W, typename System>
std::vector<tangent_row> sweep_tangent(const System& start,
                                       const OdeControl& control, Method method,
                                       std::span<const schedule_row> schedule,
                                       const tangent_directions& wrt,
                                       std::span<const std::size_t> rows);

template <typename System>
class Solver
{
public:
  using value_type = typename System::value_type;

  Solver(System sys_, OdeControl control, Method method = Method::rkck)
    : system(sys_), control_(control), solver(system, control, method)
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

  // The size of the step that reached each time in times(); NaN for the first,
  // which no step reached.
  std::vector<double> step_sizes() const { return solver.get_step_sizes(); }

  // The System this solver steps, so a caller taking the solver can constrain on
  // it rather than on what get_system_ref() happens to return.
  using system_type = System;

  System get_system() const { return system; }
  System& get_system_ref() { return system; }

  // The control this solver was built with, so a driver builds the active solver
  // with the same integration settings.
  const OdeControl& control() const { return control_; }

  // The stepper this solver was built with, fixed for its life; a recording
  // replays only under the one that made it.
  Method method() const { return solver.get_method(); }

  // Synchronize internal ODE buffers from the current system state without
  // resetting solver history/step-size state.
  void set_state_from_system()
  {
    solver.set_state_from_system(system);
  }

  // Record the wider state an insertion just reached, on the row it followed.
  void push_insertion() { solver.push_insertion(system); }

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

  // Step over a schedule, landing on each of its times: at the recorded size
  // where one is known, and to the time itself where it is not. See schedule_row
  // for why neither is derived from the other.
  //
  // A row marked an insertion is followed by the System's own state map, applied
  // here. A schedule that says where its insertions are is one a walk can execute
  // without a segment loop wrapped around it -- and it is the same map the sweep
  // transposes, so a tangent run through here traverses exactly the function
  // under test rather than a second spelling of it.
  void advance_schedule(const std::vector<ode::schedule_row>& schedule)
  {
    if (schedule.empty())
    {
      util::stop("'schedule' must hold at least the row it starts from");
    }
    if (schedule.front().insertion || !std::isnan(schedule.front().step_size))
    {
      util::stop("A schedule's first row is where it starts, which no step "
                 "reached, so it must be a step of NaN size");
    }

    // Held across the walk rather than made per insertion. `widened` is written
    // and not read: what the map leaves on the System is what this wants.
    std::vector<value_type> before;
    std::vector<value_type> widened;

    if (collect)
    {
      history.push_back(system);
    }

    for (std::size_t k = 1; k < schedule.size(); ++k)
    {
      if (schedule[k].insertion)
      {
        before.assign(system.ode_size(), value_type(0.0));
        system.ode_state(before.begin());
        ode::apply_insertion(system, schedule[k].time, before.begin(), widened);
        set_state_from_system();
        solver.push_insertion(system);
        continue;
      }
      if (std::isnan(schedule[k].step_size))
      {
        solver.step_to(system, schedule[k].time);
      }
      else
      {
        solver.step_by(system, schedule[k]);
      }
      if (collect)
      {
        history.push_back(system);
      }
    }
  }

  // The same walk over a RECORDING, which is a schedule plus what each step left.
  // Each step is taken at the size that run took, and its stages
  // LOAD what that run solved for instead of solving again -- so a pass that must
  // not re-decide (anything re-running the model to tape it) traverses the
  // function the run computed rather than a second spelling of it.
  //
  // The row travels WITH the step, because `step_record` is the schedule row plus
  // what the step left.
  //
  // ⚠️ `rec[k].solved` is what the stages of the step that REACHED `rec[k].time`
  // solved, which is the same pairing `sweep_adjoint` walks. Off by one here and
  // every stage loads its neighbour's answer, finitely.
  void replay(std::span<const ode::step_record<System>> rec)
  {
    if (rec.empty())
    {
      util::stop("'rec' must hold at least the row it starts from");
    }
    if (rec.front().insertion || !std::isnan(rec.front().step_size))
    {
      util::stop("A recording's first row is where it starts, which no step "
                 "reached, so it must be a step of NaN size");
    }

    std::vector<value_type> before;
    std::vector<value_type> widened;

    if (collect)
    {
      history.push_back(system);
    }

    for (std::size_t k = 1; k < rec.size(); ++k)
    {
      if (rec[k].insertion)
      {
        before.assign(system.ode_size(), value_type(0.0));
        system.ode_state(before.begin());
        ode::apply_insertion(system, rec[k].time, before.begin(), widened);
        set_state_from_system();
        solver.push_insertion(system);
        continue;
      }
      if (std::isnan(rec[k].step_size))
      {
        // A recording's steps are steps that were taken, so every one of them
        // has a size. A NaN here is a grid someone built by hand and called a
        // recording, and stepping TO the time would leave the rows unread.
        util::stop("A recorded step carries the size it took; entry " +
                   util::to_string(k) + " has none, so it is a grid rather "
                   "than a recording and cannot supply what its stages solved");
      }
      if (rec[k].subdivided)
      {
        util::stop("A recorded step that was subdivided cannot be replayed: "
                   "entry " + util::to_string(k) + " reached t=" +
                   util::format_double(rec[k].time) + " in several sub-steps "
                   "and holds only the last one's solved values");
      }
      if (rec[k].method != method())
      {
        util::stop("A recorded step replays only under the stepper that took "
                   "it: entry " + util::to_string(k) + " was taken by method='" +
                   method_name(rec[k].method) + "' and this solver steps by "
                   "method='" + method_name(method()) + "'. Build the solver "
                   "with the recording's method, or step steps_taken() with "
                   "advance_schedule(), which "
                   "carries no rows, to solve the stages afresh");
      }
      solver.step_by(system, rec[k], &rec[k].solved);
      if (collect)
      {
        history.push_back(system);
      }
    }
  }

  // Take a series of plain forward-Euler steps over the supplied grid. One
  // derivative evaluation per step, no error control (cf. advance_fixed, which
  // drives the selected stepper). Collects history at each supplied time.
  void advance_euler(std::vector<double> times)
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
      solver.step_euler(system, *t++);
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

  // One adaptive step that stops at time_max (Inf for no bound); see
  // SolverInternal::step(System&, double).
  void step(double time_max)
  {
    solver.step(system, time_max);
    if (collect)
    {
      history.push_back(system);
    }
  }

  // The state at each of `times`, one entry per time, the first of which must
  // be the current time. With `dense` the steps are the controller's own and
  // each requested time is read off the interpolant of the step that spans it:
  // the integration costs what it costs, however many rows are asked
  // for, and the last time is still landed on exactly. Without it every
  // requested time is landed on, as advance_adaptive() does.
  std::vector<ode::state_type<System>> advance_collect(const std::vector<double>& times,
                                                       bool dense = true)
  {
    if (times.empty())
    {
      util::stop("'times' must be vector of at least length 1");
    }
    if (!util::identical(times[0], time()))
    {
      util::stop("First element in 'times' must be same as current time");
    }
    for (size_t k = 1; k < times.size(); ++k)
    {
      if (!(times[k] > times[k - 1]))
      {
        util::stop("'times' must be strictly increasing");
      }
    }
    std::vector<ode::state_type<System>> out;
    out.reserve(times.size());
    out.push_back(state());
    if (!dense)
    {
      std::vector<double> leg(2);
      for (size_t k = 1; k < times.size(); ++k)
      {
        leg[0] = time();
        leg[1] = times[k];
        advance_adaptive(leg);
        out.push_back(state());
      }
      return out;
    }
    const double final_time = times.back();
    ode::state_type<System> y_at;
    for (size_t k = 1; k < times.size(); ++k)
    {
      const double from = time();
      for (size_t taken = 0; time() < times[k]; ++taken)
      {
        solver.check_max_steps(taken, from, times[k]);
        step(final_time);
      }
      if (util::identical(time(), times[k]))
      {
        out.push_back(state());
      }
      else
      {
        solver.interpolate(times[k], y_at);
        out.push_back(y_at);
      }
    }
    return out;
  }

  // Step adaptively until g(t, y) reaches zero or the time reaches time_max,
  // whichever comes first, and return the time g reached zero, or nothing
  // where time_max came first.
  //
  // An event is a change in the sign g had where the call started, seen at
  // the end of a step. It is located on that step's interpolant (locate(), no
  // rates evaluated) to `tol` in time, the step is taken back
  // (retract_step()), and the solver then steps to the located time. The state
  // left there is therefore an integrated one, to the solver's own tolerance
  // under every stepper. It and everything after it are what a run would
  // hold that took the same steps up to the one that crossed and then asked to
  // land on the event. g = 0 at the start is an event at the start, and
  // stepping stops there.
  //
  // A step whose ends share a sign is taken to cross nothing, so a g that
  // crosses and crosses back inside one step goes unseen; step_size_max is
  // the bound on how short-lived a crossing can be and still be found.
  //
  // At a located event g is zero only to the location's tolerance and the
  // solver's, so a caller that asks again with the same g starts from a value
  // of either sign near zero, and can be stopped again straight away. Ask for
  // the next level instead, as a sequence of sizes does.
  //
  // `g` is called as g(t, y) with y a const state_type&, and must not change
  // the System: the interpolant's states are not ones the System holds. A
  // non-finite g is refused. Anything thrown while the crossing is located
  // takes the crossing step back first, so the solver is left between steps,
  // at the start of that step.
  template <class G>
  std::optional<double> advance_until(G&& g, double time_max,
                                      const root::tolerance& tol)
  {
    if (!(std::isfinite(time_max) && time_max >= time()))
    {
      util::stop("advance_until(): time_max = " + util::format_double(time_max) +
                 " must be finite and not before the current time " +
                 util::format_double(time()));
    }
    const auto g_finite = [&g](double t, const ode::state_type<System>& y) {
      const double value = g(t, y);
      if (!std::isfinite(value))
      {
        util::stop("advance_until(): g is not finite at t = " +
                   util::format_double(t));
      }
      return value;
    };
    double g_prev = g_finite(time(), solver.state_ref());
    if (g_prev == 0.0)
    {
      return time();
    }
    const bool positive = g_prev > 0.0;
    const double from = time();
    for (size_t taken = 0; time() < time_max; ++taken)
    {
      solver.check_max_steps(taken, from, time_max);
      solver.hold_step_start();
      solver.step(system, time_max);
      const double g_end = g_finite(time(), solver.state_ref());
      if (g_end == 0.0)
      {
        if (collect)
        {
          history.push_back(system);
        }
        return time();
      }
      if ((g_end > 0.0) == positive)
      {
        g_prev = g_end;
        if (collect)
        {
          history.push_back(system);
        }
        continue;
      }
      double t_event;
      try
      {
        t_event = solver.locate(g_finite, g_prev, g_end, tol);
      }
      catch (...)
      {
        solver.retract_step(system);
        throw;
      }
      solver.retract_step(system);
      while (time() < t_event)
      {
        solver.check_max_steps(taken++, from, time_max);
        step(t_event);
      }
      return t_event;
    }
    return std::nullopt;
  }

  double get_step_size() const { return solver.get_step_size(); }
  void set_step_size(double h) { solver.set_step_size(h); }
  std::size_t get_n_rejections() const { return solver.get_n_rejections(); }
  bool mid_step() const { return solver.mid_step(); }
  const SolverInternal<System>& get_internal() const { return solver; }

  bool get_collect() const { return collect; }

  void set_collect(bool x) { collect = x; }

  std::size_t get_history_size() const { return history.size(); }

  std::vector<System> get_history() const { return history; }

  System get_history_step(std::size_t i) const { return history.at(i); }

  // Keep the state at each accepted step beside the time and the size that reached
  // it. Set before the first step; the row the run starts from is filled from the
  // state the solver holds when this is turned on. set_state() and reset() start
  // a new record.
  void set_keep_states(bool keep) { solver.set_keep_states(keep); }
  // The record the run kept: one row per accepted step, each carrying its time,
  // the size that reached it and the state there. What a sweep reads, and the
  // only place it reads them from.
  std::span<const ode::step_record<System>> recording() const {
    return solver.recording();
  }
  // The steps this run took, without its insertion rows. Read off the record, so
  // the time and the size that reached it cannot be paired across two runs.
  std::vector<ode::schedule_row> steps_taken() const {
    return solver.steps_taken();
  }

  // Carry lambda back over recorded rows k_last down to k_first + 1, highest
  // first, adding each row's parameter contribution into `parameter_adjoint`.
  //
  // Wherever the state changed width inside those rows the sweep stops, narrows
  // the System across it and transposes the map that widened it, so the rows it
  // carries arrive by the same route a step's do. A System that declares
  // assign_from serves the whole sweep as one active System: rebound once, and
  // re-seated at each width after that (adjoint_system::reseat). Any other is
  // rebound at each width, as a fresh active System.
  //
  // An output adjoint that is a multiple of an earlier one is carried down as
  // that one and rebuilt at the end (adjoint_multiples), so it costs no sweep of
  // its own. Its state and parameter adjoints then agree with a sweep of its own
  // to rounding; every other output's are what they would be swept alone.
  //
  // `extra_stops` names further rows to stop at. The adjoint carried across such a
  // stop is the same one either way, so a caller can ask for a split and compare: a
  // check, not a choice.
  //
  // The System is put where each segment needs it and left where the run left it,
  // so a sub-segment can be swept without the caller positioning anything and the
  // whole can be swept again afterwards.
  //
  // Returns how many segments carried a step.
  std::size_t sweep_adjoint(ode::adjoint_batch& lambda,
                            ode::adjoint_batch& parameter_adjoint,
                            size_t k_first, size_t k_last,
                            const std::vector<size_t>& extra_stops = {})
  {
    if (&lambda == &parameter_adjoint) {
      util::stop("sweep_adjoint: the state adjoints are replaced and the "
                 "parameter adjoints accumulated, so they cannot be the same "
                 "batch");
    }
    const ode::adjoint_multiples multiples = ode::adjoint_multiples::of(lambda);
    if (!multiples.any()) {
      return sweep_rows(lambda, parameter_adjoint, k_first, k_last, extra_stops);
    }
    // Asked here, before the batch is narrowed: the recordings below see only
    // the kept rows, and the rebuild writes every row.
    if (parameter_adjoint.size() < lambda.size()) {
      util::stop("sweep_adjoint: one row of parameter adjoints per output, to "
                 "accumulate into");
    }
    // The kept rows' parameter adjoints start where the caller's do and are
    // written back whole, so they accumulate in the order a sweep of the whole
    // batch would; a multiple takes c times what its row added.
    ode::adjoint_batch carried = lambda.select(multiples.kept);
    const ode::adjoint_batch before = parameter_adjoint.select(multiples.kept);
    ode::adjoint_batch carried_parameters = before;
    const std::size_t swept = sweep_rows(carried, carried_parameters, k_first,
                                         k_last, extra_stops);
    lambda.assign(lambda.size(), carried.width());
    for (std::size_t m = 0; m < lambda.size(); ++m) {
      const std::size_t b = multiples.from[m];
      const double c = multiples.factor[m];
      const std::span<double> state = lambda[m], parameters = parameter_adjoint[m];
      if (multiples.kept[b] == m) {
        std::copy(carried[b].begin(), carried[b].end(), state.begin());
        std::copy(carried_parameters[b].begin(), carried_parameters[b].end(),
                  parameters.begin());
        continue;
      }
      for (std::size_t i = 0; i < state.size(); ++i) {
        state[i] = c * carried[b][i];
      }
      for (std::size_t p = 0; p < parameters.size(); ++p) {
        parameters[p] += c * (carried_parameters[b][p] - before[b][p]);
      }
    }
    return swept;
  }

  // The whole recording.
  std::size_t sweep_adjoint(ode::adjoint_batch& lambda,
                            ode::adjoint_batch& parameter_adjoint,
                            const std::vector<size_t>& extra_stops = {})
  {
    return sweep_adjoint(lambda, parameter_adjoint, 0, recording().size() - 1,
                         extra_stops);
  }

  // The forward counterpart of sweep_adjoint: this run's recording walked again
  // from its first row, carrying one direction for each input `wrt` names, W
  // at a time. Returns the state and its derivative along each direction at
  // each of `rows`, indices into recording(), every row when none are given.
  // The steps are the run's, so over the same rows this and a sweep
  // differentiate the same function and agree to rounding.
  //
  // A subdivided row is refused, as a sweep refuses it: its size is the
  // interval's, and one step across the interval is not what the run took.
  template <std::size_t W>
  std::vector<ode::tangent_row> sweep_tangent(const ode::tangent_directions& wrt,
                                              std::vector<std::size_t> rows = {})
      const
  {
    const auto rec = recording();
    if (rec.empty())
    {
      util::stop("sweep_tangent: nothing recorded; set_keep_states(true) "
                 "before the run");
    }
    for (std::size_t k = 0; k < rec.size(); ++k)
    {
      if (rec[k].subdivided)
      {
        util::stop("sweep_tangent: recorded step " + util::to_string(k) +
                   " was subdivided, so its row is not one step of the run");
      }
    }
    if (rows.empty())
    {
      rows.resize(rec.size());
      for (std::size_t k = 0; k < rows.size(); ++k)
      {
        rows[k] = k;
      }
    }
    // Only as far as the last row asked for: nothing after it is read.
    const std::size_t end = *std::max_element(rows.begin(), rows.end()) + 1;
    if (end > rec.size())
    {
      util::stop("sweep_tangent: row " + util::to_string(end - 1) +
                 " is past the recording's " + util::to_string(rec.size()));
    }
    const std::vector<ode::schedule_row> schedule(rec.begin(), rec.begin() + end);
    System start = system;
    ode::be_at_step(start, rec, 0);
    return ode::sweep_tangent<W>(start, control_, method(), schedule, wrt, rows);
  }

  // Rate evaluations the sweeps since the last clear have recorded;
  // Step::recorded_rates says what the number means and what it catches.
  std::size_t recorded_rates() const { return solver.recorded_rates(); }
  void clear_recorded_rates() { solver.clear_recorded_rates(); }

  // Should we record history at every step?
  // TODO: should this be part of ode_solver?
std::vector<System> history;

private:
  // sweep_adjoint, every output adjoint carried down.
  std::size_t sweep_rows(ode::adjoint_batch& lambda,
                         ode::adjoint_batch& parameter_adjoint,
                         size_t k_first, size_t k_last,
                         const std::vector<size_t>& extra_stops)
  {
    using scalar = ode::adjoint_scalar<double>;
    static_assert(ode::Sweepable<System>,
                  "this System cannot be swept: it needs rebind<U>(), "
                  "ad_parameters(), for_each_active(f) and "
                  "set_recorded_state(y, time); see Sweepable in "
                  "ode_interface.hpp");
    const std::span<const ode::step_record<System>> rec = recording();
    if (rec.size() < 2) {
      util::stop("sweep_adjoint: no recorded steps to sweep; run the adaptive "
                 "pass first");
    }
    if (k_first >= k_last || k_last >= rec.size()) {
      util::stop("sweep_adjoint: the rows to sweep are not recorded rows");
    }

    // A insertion is a row the sweep carries an adjoint ACROSS, so it cannot be one
    // the sweep starts at or the one it is left standing on. Neither happens on a
    // recording a run made -- the schedule puts every insertion at the start of
    // an interval that then steps -- and both are checked here rather than at each
    // be_at_step, because the width this leaves the System at is a promise and the
    // call that restores it cannot raise.
    if (rec.back().insertion || rec[k_last].insertion) {
      util::stop("sweep_adjoint: a recording cannot end at an insertion, because "
                 "nothing stepped away from the state it made");
    }

    // The width on exit is a promise, and a throw is an exit, so a destructor
    // keeps it. The sweep starts at the run's own width and narrows as it goes,
    // so a sweep abandoned high up leaves the System at its widest -- where every
    // caller's tail widens back from the lowest and reads the mismatch as a length
    // error one call later, naming neither this walk nor what refused.
    //
    // The solver's own buffers are reset from the System at the same time:
    // a sweep sizes the stage buffers to each segment it walks, so after a segment
    // narrower than the run they hold a truncated state marked current, and a
    // step taken next would start from it.
    struct restore_on_exit {
      System& sys;
      SolverInternal<System>& solver;
      std::span<const ode::step_record<System>> rec;
      ~restore_on_exit() {
        // This runs with another exception possibly in flight, so a failure here
        // cannot be raised: it would end the process rather than the call that is
        // already failing.
        try {
          ode::be_at_step(sys, rec, rec.size() - 1);
          solver.set_state_from_system(sys);
        } catch (...) {
        }
      }
    } restore{system, solver, rec};

    // One tape for the whole sweep, held active across every recording it takes.
    // Clearing between recordings keeps the capacity the largest of them grew,
    // where one tape per recording regrows it every time.
    ode::adjoint_tape<double> tape(false);
    ode::tape_scope<ode::adjoint_tape<double>> running{tape};

    // The active System every recording below is taken on. A System that
    // declares assign_from is built at the first and re-seated at each after
    // it, so a sweep pays one System construction however many widths it
    // passes through; any other is built afresh at each. Declared after the
    // tape's scope, so it releases its slots with the tape still active.
    std::optional<ode::adjoint_system<System>> active;
    const auto seat = [&]() -> ode::adjoint_system<System>& {
      if (active) {
        if constexpr (ode::AssignsFrom<
                          typename ode::adjoint_system<System>::system_type,
                          System>) {
          active->reseat(system);
          return *active;
        } else {
          active.reset();
        }
      }
      active.emplace(system, tape);
      return *active;
    };

    // One segment per width, highest first: the System is put at the segment's
    // top and the active System seated on it there, so no sweep inherits a
    // width another call left behind.
    //
    // A stop is either an insertion the run recorded or a cut a caller asked for, and
    // they leave the sweep in different places. A cut is a row the sweep resumes
    // at; an insertion is a row the sweep carries the adjoint ACROSS, so it resumes
    // one row below, on the state the map ran on.
    size_t swept = 0;
    size_t hi = k_last;
    const auto sweep_down_to = [&](size_t lo) -> void {
      ode::be_at_step(system, rec, hi);
      // A segment with no step in it cuts nothing, which is what two stops in a row
      // gives. Nor does a segment with no state -- a run that starts empty and
      // widens at a later insertion, as plant's invader does before its first
      // cohort: a step with no state has no output, so it carries no adjoint and
      // adds nothing to any parameter's.
      if (lo < hi && lambda.width() > 0) {
        sweep_segment(seat(), rec, lambda, parameter_adjoint, lo, hi);
        ++swept;
      }
    };

    // The rows this can stop at are the rows strictly inside the segment, which is
    // what the bounds say, and each is reached once. A row that is both an insertion
    // and a cut is the insertion, which is what keeps a cut free of a map.
    for (size_t at = k_last; at-- > k_first + 1;) {
      const bool insertion = rec[at].insertion;
      if (!insertion && std::find(extra_stops.begin(), extra_stops.end(), at) ==
                           extra_stops.end()) {
        continue;
      }
      sweep_down_to(at);
      if (!insertion) {
        hi = at;
        continue;
      }
      // An insertion that leaves the state empty -- one that introduced nothing
      // into a run with no state yet, as a replay's earlier introductions do for
      // an invader whose own come later -- has no output to transpose, and the
      // state below it is empty too.
      if (lambda.width() == 0) {
        hi = at - 1;
        continue;
      }
      // The map that widened `at`, transposed at the width below it -- which is
      // the width the segment below runs at, and the width the row below holds.
      //
      // Applying the map widens the active System it is recorded on, so the
      // segment below re-seats it at its own width before recording anything.
      ode::be_at_step(system, rec, at - 1);
      const double when = rec[at].time;
      // A System that declares no map passes its state through, which is only a
      // map where the width did not move. Checked here, because the pass-through
      // reads one entry per output and the input it would read past is the
      // narrower state.
      if constexpr (!requires(System& s, double t,
                              typename state_type<System>::const_iterator x,
                              state_type<System>& y) {
                      s.apply_insertion(t, x, y);
                    }) {
        if (rec[at].state.size() != rec[at - 1].state.size()) {
          util::stop("sweep_adjoint: the recording widens from " +
                     util::to_string(static_cast<int>(rec[at - 1].state.size())) +
                     " to " +
                     util::to_string(static_cast<int>(rec[at].state.size())) +
                     " at row " + util::to_string(at) +
                     ", but the System declares no apply_insertion, so there "
                     "is no map to transpose");
        }
      }
      auto insert = [&](auto& sys,
                        typename std::vector<scalar>::const_iterator x,
                        std::vector<scalar>& y) -> void {
        ode::apply_insertion(sys, when, x, y);
      };
      ode::adjoint_system<System>& widened = seat();
      ode::adjoint_batch narrowed;
      ode::state_and_parameter_adjoints(widened, rec[at - 1].state, lambda, insert,
                                       narrowed, parameter_adjoint);
      lambda = std::move(narrowed);
      hi = at - 1;
    }
    sweep_down_to(k_first);
    return swept;
  }

  // One segment, on the active System the caller seated at this segment's width:
  // every step in it is recorded on that one System, each recording releasing its
  // slots before taking the next.
  //
  // Every state visited has to be that width, so a segment an insertion crossed is
  // refused by the length check inside the loop rather than swept at one width.
  void sweep_segment(ode::adjoint_system<System>& active,
                     std::span<const ode::step_record<System>> rec,
                     ode::adjoint_batch& lambda,
                     ode::adjoint_batch& parameter_adjoint,
                     size_t k_first, size_t k_last)
  {
    if (lambda.empty()) {
      util::stop("sweep_adjoint: needs at least one output adjoint");
    }
    // Asked of the System the recordings are taken on, which is the width every
    // state below has to be loaded at. Named, because a bare length mismatch
    // here is read as the caller's and says nothing about the seam it is really
    // about: a batch carried at one width against a System at another. One
    // width for every output adjoint, so this is asked of the batch and not of
    // each adjoint in it.
    if (lambda.width() != active.system.ode_size()) {
      util::stop("sweep_adjoint: the output adjoints are " +
                 util::to_string(static_cast<int>(lambda.width())) +
                 " wide against a System of " +
                 util::to_string(static_cast<int>(active.system.ode_size())) +
                 ", so the two are not at the same insertion");
    }
    // ⚠️ THE DESCENT STOPS AT THE FIRST NON-FINITE ENTRY, and must not be changed
    // to carry it. This hands its caller one number per input, so an overflow the
    // sweep picked up thousands of steps earlier would arrive indistinguishable
    // from a NaN the last step made -- and a caller polling for a declared
    // refusal would find none, which is the one failure this gradient must never
    // produce. Raised as AdjointOverflowError, so a caller can tell it from any
    // other failure.
    //
    // ODELIA_ADJOINT_TRACE=steps prints the magnitude at every step, which is
    // what says whether the sweep compounded into the failure or met it. One
    // line per recorded step is thousands of them, so it is not the default; the
    // magnitude one step above is in the message either way.
    static const char* const trace = std::getenv("ODELIA_ADJOINT_TRACE");
    double worst_before = 0.0;
    const bool per_step = trace != nullptr && std::strcmp(trace, "steps") == 0;
    ode::adjoint_batch lambda_in;
    for (size_t k = k_last; k > k_first; --k) {
      // What the run's step k ran from: the row below it, whether that row is a
      // step's landing or an insertion's output.
      const state_type<System>& from = rec[k - 1].state;
      util::check_length(from.size(), active.system.ode_size());
      if (rec[k].subdivided) {
        util::stop("sweep_adjoint: recorded step " + util::to_string(k) +
                   " reached t=" + util::format_double(rec[k].time) +
                   " in several sub-steps, and one step of its interval is "
                   "not the step the run took; record the run adaptively, or "
                   "pin it finely enough that nothing is refused");
      }
      // The run's k1, taken as step_record::k1_time says the run took it. No row
      // holds a k1 derived afresh, so a System that stands in a loaded value
      // across rows has to supply it there itself.
      const bool k1_carried = !std::isnan(rec[k].k1_time);
      if (k1_carried && (k < 2 || rec[k - 1].insertion)) {
        util::stop("sweep_adjoint: recorded step " + util::to_string(k) +
                   " says its first rates were carried from the step below, "
                   "but row " + util::to_string(k - 1) + " is not a step");
      }
      const ode::carried_rates<System> k1{rec[k - 1].solved, rec[k].k1_time};
      solver.step_adjoint(active, rec[k].solved, k1_carried ? &k1 : nullptr,
                          rec[k - 1].time, rec[k].step_size, from,
                          lambda, lambda_in, parameter_adjoint,
                          rec[k].subsystem ? &*rec[k].subsystem : nullptr);
      // Swapped rather than moved from: a move leaves the buffer this step wrote
      // into empty, so the next step allocates one the same size again. Swapping
      // hands it the row above's, which the sweep refills rather than regrows.
      std::swap(lambda, lambda_in);
      double worst_here = 0.0;
      for (size_t m = 0; m < lambda.size(); ++m) {
        const std::span<double> row = lambda[m];
        for (size_t j = 0; j < row.size(); ++j) {
          if (!std::isfinite(row[j])) {
            util::stop_adjoint_overflow(
                "the adjoint left the representable range at step " +
                util::to_string(k) + " of [" + util::to_string(k_first) + ", " +
                util::to_string(k_last) + "], t=" +
                util::format_double(rec[k - 1].time) + ", h=" +
                util::format_double(rec[k].step_size) + ": output adjoint " +
                util::to_string(m) + "'s entry " + util::to_string(j) + " of " +
                util::to_string(row.size()) + " is " +
                util::format_double(row[j]) + ", where the step above carried " +
                util::format_double(worst_before) +
                ". A sweep is a product of step Jacobians and has no error "
                "control, so it can pass outside the range of the answer it "
                "returns; this one did not come back.");
          }
          const double at = std::fabs(row[j]);
          if (at > worst_here) worst_here = at;
        }
      }
      if (per_step) {
        std::fprintf(stderr,
                     "ODELIA_ADJOINT_TRACE step %zu of [%zu, %zu] t=%.12g "
                     "h=%.6g worst |lambda| %.6g\n",
                     k, k_first, k_last, rec[k - 1].time, rec[k].step_size,
                     worst_here);
      }
      worst_before = worst_here;
    }
  }

  bool collect;
  System system;
  OdeControl control_;
  SolverInternal<System> solver;
};

// The forward counterpart of a sweep: `schedule` walked from `start`, which
// stands where the schedule starts, with the System rebound to a scalar
// carrying one lane per direction `wrt` names. Returns the state and its
// derivative along each direction at each of `rows`, indices into the schedule,
// in the order asked.
//
// The steps are the schedule's, at its sizes, and nothing a lane carries reaches
// a step's size or whether it is taken, so the derivative is that of the scheme
// as those steps took it -- the convention sweep_adjoint keeps. Under RODAS the
// Jacobian carries the lanes too (Jacobian::compute), so a step's derivative is
// the step's own. A row of NaN size is stepped to its time and may be
// subdivided there, which is differentiated as taken.
//
// More directions than W are walked in further passes of W. The lanes do not
// mix, so the answer does not depend on W; the value is the same in every pass.
template <std::size_t W, typename System>
std::vector<tangent_row> sweep_tangent(const System& start,
                                       const OdeControl& control, Method method,
                                       std::span<const schedule_row> schedule,
                                       const tangent_directions& wrt,
                                       std::span<const std::size_t> rows)
{
  using lane_type = tangent_lanes<W>;
  static_assert(W > 0, "sweep_tangent needs at least one lane");
  static_assert(TangentSweepable<System, lane_type>,
                "sweep_tangent needs a System that rebinds to the lane scalar, "
                "stands on a given state there (set_recorded_state), and "
                "solves for no inner values (see TangentSweepable)");
  using lane_system = typename rebound_system<System, lane_type>::type;

  if (schedule.empty())
  {
    util::stop("sweep_tangent: 'schedule' must hold at least the row it "
               "starts from");
  }
  if (method == Method::rodas && !RodasStep<lane_system>::supported)
  {
    util::stop("sweep_tangent: under method='rodas' the System's Jacobian is a "
               "tangent of a tangent, which it must declare it compiles "
               "(`static constexpr bool nests_tangents = true`; see "
               "NestsTangents) unless it has its own ode_jacobian()");
  }
  for (const std::size_t r : rows)
  {
    if (r >= schedule.size())
    {
      util::stop("sweep_tangent: row " + util::to_string(r) +
                 " is past the schedule's " + util::to_string(schedule.size()));
    }
  }
  const std::size_t n_start = start.ode_size();
  for (const std::size_t i : wrt.initial)
  {
    if (i >= n_start)
    {
      util::stop("sweep_tangent: initial-state direction " + util::to_string(i) +
                 " is past the state's " + util::to_string(n_start) +
                 " components");
    }
  }
  std::vector<double> y_start(n_start);
  start.ode_state(y_start.begin());
  const double t_start = schedule.front().time;

  const std::vector<schedule_row> walk(schedule.begin(), schedule.end());
  const std::size_t n_parameters = wrt.parameters.size();
  const std::size_t n_directions = wrt.size();
  std::vector<tangent_row> out(rows.size());
  for (tangent_row& row : out)
  {
    row.directions = n_directions;
  }

  for (std::size_t first = 0; first == 0 || first < n_directions; first += W)
  {
    const std::size_t last = std::min(first + W, n_directions);
    lane_system walker = start.template rebind<lane_type>();

    // Parameters before the state: the state's load computes the rates the
    // first step starts from, and they must carry the parameters' lanes.
    if (first < n_parameters)
    {
      if constexpr (has_ad_parameters<lane_system>::value)
      {
        const std::vector<lane_type*> parameters = walker.ad_parameters();
        for (std::size_t k = first; k < std::min(last, n_parameters); ++k)
        {
          const std::size_t p = wrt.parameters[k];
          if (p >= parameters.size())
          {
            util::stop("sweep_tangent: parameter direction " +
                       util::to_string(p) + " is past the System's " +
                       util::to_string(parameters.size()) + " ad_parameters()");
          }
          set_direction(*parameters[p], k - first, 1.0);
        }
      }
      else
      {
        util::stop("sweep_tangent: a direction on a parameter needs "
                   "ad_parameters() on the System");
      }
    }
    std::vector<lane_type> y(n_start);
    for (std::size_t i = 0; i < n_start; ++i)
    {
      y[i] = lane_type(y_start[i]);
    }
    for (std::size_t k = std::max(first, n_parameters); k < last; ++k)
    {
      set_direction(y[wrt.initial[k - n_parameters]], k - first, 1.0);
    }
    walker.set_recorded_state(y, t_start);

    Solver<lane_system> s(walker, control, method);
    s.set_collect(false);
    s.set_keep_states(true);
    s.advance_schedule(walk);
    const auto rec = s.recording();
    if (rec.size() != schedule.size())
    {
      util::stop("sweep_tangent: the walk kept " + util::to_string(rec.size()) +
                 " rows of a schedule of " + util::to_string(schedule.size()));
    }

    for (std::size_t j = 0; j < rows.size(); ++j)
    {
      const auto& state = rec[rows[j]].state;
      tangent_row& row = out[j];
      if (first == 0)
      {
        row.value.resize(state.size());
        row.derivative.assign(state.size() * n_directions, 0.0);
        for (std::size_t i = 0; i < state.size(); ++i)
        {
          row.value[i] = xad::value(state[i]);
        }
      }
      for (std::size_t i = 0; i < state.size(); ++i)
      {
        for (std::size_t k = first; k < last; ++k)
        {
          row.derivative[i * n_directions + k] =
              derivative_along(state[i], k - first);
        }
      }
    }
  }
  return out;
}
}
}
#endif
