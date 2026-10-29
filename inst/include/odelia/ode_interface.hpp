// -*-c++-*-
#ifndef ODELIA_ODE_INTERFACE_HPP_
#define ODELIA_ODE_INTERFACE_HPP_

#include <cstddef>
#include <algorithm>
#include <span>
#include <array>
#include <optional>
#include <vector>
#include <odelia/adjoint_scalar.hpp>
#include <odelia/ode_util.hpp>
#include <odelia/tangent.hpp>

#include <concepts>
#include <iterator>
#include <limits>
#include <type_traits>
#include <utility>
#include <array>
#include <vector>
#include <XAD/XAD.hpp>

namespace odelia {
namespace ode {

// Aliases for callers written against the double-only state. The templated
// state_type below is the one the solver uses.
typedef std::vector<double>::const_iterator const_iterator;
typedef std::vector<double>::iterator       iterator;

// Type alias for state vectors based on System's value_type
template<typename System>
using state_type = std::vector<typename System::value_type>;

// A System that carries its own clock. One that does not is time homogeneous,
// and the calls below hand it no time rather than refusing it.
template <typename T>
concept HasOdeTime = requires(const T& t) {
  { t.ode_time() } -> std::convertible_to<double>;
};

// A System that can hand back a copy of itself on scalar U. The rebound type must
// itself be a System on U: a rebind() returning the wrong scalar fails here
// rather than on the first arithmetic inside the caller.
//
// U is named rather than defaulted. Defaulting it to the System's own scalar asks
// whether a System can rebind to the scalar it already has, which is a different
// question from the one every caller means and is answered yes by types that
// cannot do what the caller needs.
template <typename S, typename U>
concept Rebindable = requires(const S& s) {
  { s.template rebind<U>() };
  requires std::same_as<typename decltype(s.template rebind<U>())::value_type, U>;
};

// The System type rebound to scalar U, i.e. decltype(system.rebind<U>()).
// When the System has no rebind() the type is not evaluated (a harmless
// placeholder is used instead), so a holder can still be class-instantiated for
// systems that will never rebind; the actual use is gated on the concept.
template <typename S, typename U, bool = Rebindable<S, U>>
struct rebound_system {
  using type = decltype(std::declval<const S>().template rebind<U>());
};
template <typename S, typename U>
struct rebound_system<S, U, false> {
  using type = S;
};

// What one rate evaluation solved for that its state does not determine: the
// results of the inner solves it ran, in the order it ran them. A System that
// solves for nothing declares nothing and this is empty.
//
// Only values of that kind belong here. Anything a rate evaluation can recompute
// from the state it was handed MUST be recomputed, because the transpose runs
// through it -- loading such a value as a recorded double severs the tape there
// and the sweep comes back wrong with every number finite.
struct no_solved_values {};

template <class System>
struct solved_values {
  using type = no_solved_values;
};
template <class System>
  requires requires { typename System::solved_values; }
struct solved_values<System> {
  using type = typename System::solved_values;
};
template <class System>
using solved_values_t = typename solved_values<System>::type;

// A System whose rate evaluation solves for something its state does not
// determine: the branch of an inner root-find, an early exit, the point an
// optimisation landed on. A pass re-running the model in order to tape it has to
// take the RUN'S result rather than solve again -- re-solving risks landing the
// other side of a branch, and what is then taped is a function the run never
// computed, with every number finite.
//
// So the run STORES what it solved for and a later pass LOADS it. Which of the two
// is happening is not a flag anyone keeps: a walk hands over the values, and
// whether it hands them over to be written or to be read is the constness of what
// it hands over.
//
// The extent is one rate evaluation, which is what `derivs` is. A walk that hands
// over nothing opens no extent, so a reload out of band cannot read a record and a
// record cannot complete a state it was not taken at.
//
// `stage` is the evaluation's place in its step's row (see solved_row), by which a
// System can find the same evaluation in another run's recording.
template <typename System>
concept SolvesForValues =
  requires(System s, solved_values_t<System>& into,
           const solved_values_t<System>& from, std::size_t stage) {
    s.store_solved(into, stage);
    s.load_solved(from, stage);
    s.end_solved();
  };

// How a step integrated the System's subsystem in substeps: the slope its inputs
// were extrapolated along, and where each substep ends, as fractions of the step.
struct subsystem_substeps {
  std::vector<double> input_slope;
  std::vector<double> ends;
  bool operator==(const subsystem_substeps&) const = default;
};

// One row of a schedule: what carries the state from one boundary to the
// next. A step reaches `time`, by `step_size` where a run pinned it and by whatever
// the controller chooses where that is NaN. A insertion applies the System's own
// state map and reaches the time it started at, taking none.
//
// Time and size are both kept because neither can be recovered from the other. A
// size differenced back out of two times is not the size that was taken, since
// fl(fl(t + h) - t) is not h; and a time reached by adding sizes is not where the
// run landed, since a run sets its last step into an interval to the interval's
// end and fl(t + (t1 - t)) is not t1. They are one object because two vectors side
// by side can be paired across different runs; one cannot.
//
// A NaN size is a time with no size known, which is a grid point rather than a
// recorded step: step TO it. So one type covers a grid a caller chose and a
// schedule a run emitted, and each row says which it is rather than a second
// container's emptiness saying it for all of them.
//
// The first row is where a schedule starts, which no step reached: a step of NaN
// size. `insertion` is last and defaults, so a grid written {time, NaN} is a
// schedule of steps without saying so.
//
// `subsystem` holds a step's substeps where it substepped the System's subsystem
// (see HasSubsystem); a pinned run and the sweep take them as recorded.
struct schedule_row {
  double time;
  double step_size;
  bool insertion = false;
  std::optional<subsystem_substeps> subsystem{};
};

// Integration method: the explicit Cash-Karp RKCK 4(5) stepper (default), the
// implicit RODAS4(3) Rosenbrock stepper for stiff systems, or the explicit
// Dormand-Prince 5(4) stepper with order-4 dense output. Here rather than with
// the solver because a recorded row names the stepper that made it.
enum class Method { rkck, rodas, dopri };

// The name R passes for each, for messages.
inline const char* method_name(Method m) {
  switch (m) {
    case Method::rodas: return "rodas";
    case Method::dopri: return "dopri";
    default: return "rkck";
  }
}

// What one step's rate evaluations solve for, in the order the explicit
// steppers (Cash-Karp and Dormand-Prince) take them. One name, because the
// forward walk, the sweep and the record (step_record::solved) all have to agree
// on the shape.
//
// Indexed by stage, counted from the second: stage s is the stepper's stage s + 2.
// SIX AND NOT FIVE. Five are the stages after the first; the sixth, `end_stage`, is
// the evaluation at the state the step ends at, which first-same-as-last hands the
// next step as its own k1. The first stage has no entry: it is carried from the
// step before or derived afresh, and step_record::k1_time says which. A SWEEP
// reads 0..4 for its stages and, where the run carried k1 over, the sixth entry
// of the step before (see carried_rates) -- and a FORWARD replay cannot
// re-derive it at all, because re-deriving is exactly what it is replaying to
// avoid, and a step whose k1 was re-derived is wrong at first order in h.
inline constexpr std::size_t end_stage = 5;
template <class System>
using solved_row = std::array<solved_values_t<System>, end_stage + 1>;

// A step's first rates where first-same-as-last carried them from the step
// before: that step's sixth evaluation, the row it was taken in and the time it
// ran at (step_record::k1_time says why that time is its own).
template <class System>
struct carried_rates {
  const solved_row<System>& row;
  double time;
};

// One row of a recording: the schedule row the run executed and the state it left
// the System holding. A recording row IS a schedule row plus its state, so it
// derives rather than repeating the two fields, and a time cannot be paired with a
// state from another container.
//
// A insertion is a row like any other, and its state is what the map produced. So
// no row carries two states and nothing has to choose between them. A insertion row
// shares its time with the row below it, which is the row that holds the state the
// map ran on.
//
// A run of n steps on a state of width w holds n rows of w scalars and one row of
// solved values each, and nothing else: stage rates and stage states are
// recomputed by the sweep, which is cheaper than holding them.
template <typename System>
struct step_record : schedule_row {
  state_type<System> state;

  // What this step's rate evaluations solved for; see solved_row.
  solved_row<System> solved;

  // The stepper that took this step. Both explicit steppers fill a row of the
  // same shape, but entry i is the answer at THEIR stage i + 2, and the two
  // tableaux put stages 4..6 at different times and states. A replay under the
  // other stepper would load every one of those at the wrong stage, finitely, so
  // replay refuses it; a schedule (steps_taken()) carries no rows and
  // runs under either.
  Method method = Method::rkck;

  // A pinned step that was refused and crossed its interval in several sub-steps.
  // The row still holds the interval, which is what a run on the caller's times
  // needs, but its solved values are the LAST sub-step's and one Runge-Kutta
  // step of the interval is not the step the run took. A sweep or a replay of
  // such a row would be wrong at every entry with every number finite, so both
  // refuse it.
  bool subdivided = false;

  // Where this step's first rates came from. First-same-as-last hands a step the
  // sixth evaluation of the step before, taken in that row's `solved[5]` at the
  // time it ran at, and that time is kept here: fl(t + h) is not always the time
  // the row below records, which a pinned interval or a replay sets exactly. NaN
  // where the run derived the rates afresh instead -- at the run's start, after
  // an insertion, after any reset from the System -- with no row at all.
  // ⚠️ A SWEEP TAKES k1 THE WAY THIS SAYS THE RUN DID, never as the row below
  // implies: a step row below does not mean the rates were carried, because a
  // reset between two steps keeps no row of its own.
  double k1_time = std::numeric_limits<double>::quiet_NaN();
};

// A System whose state vector gains entries during a run does not declare a
// concept for it. Two members carry the whole of it, and they are as mandatory as
// ode_size() for a System a sweep is asked to walk, so they are called directly
// like it:
//
//   set_recorded_state(y, time)  -- be the shape this recorded time implies, then
//                                   take these values. A run loads into the shape
//                                   it already built; a replay does not know it,
//                                   and derives it from the schedule the run was
//                                   driven by.
//   apply_insertion(time, x, y)  -- the insertion as a map, the state below it in
//                                   and the whole wider state out, so it runs at
//                                   any scalar and the sweep can transpose it. It
//                                   leaves the System holding that wider state.
//
// For a width that never moves, the first is a System's ordinary load. The second
// is asked for only where the width changed, so a System that never widens is
// never asked for it -- ode::apply_insertion below is what a walk calls, and
// passing the state through is what an insertion is for such a System.
//
// An insertion whose TIME depends on the parameters is a different map: its
// adjoint carries a term through that time which nothing here computes. A System
// walked by the sweep asserts its insertions are scheduled, not triggered.

// What a System needs beyond solving for a reverse sweep to walk it: a copy of
// itself on the adjoint scalar, the parameters the sweep accumulates adjoints
// for, every member carrying the scalar (so the tape's slots can be handed back
// before each recording is cleared), and a way to stand on a recorded state. A
// System that widens adds apply_insertion, which cannot be a concept because its
// absence is what a fixed-width System means. Solver::sweep_adjoint asserts
// this, so a System missing a member fails to compile naming the four.
template <typename S>
concept Sweepable =
  Rebindable<S, adjoint_scalar<double>> &&
  requires(typename rebound_system<S, adjoint_scalar<double>>::type a,
           const std::vector<adjoint_scalar<double>>& y, double t) {
    { a.ad_parameters() } -> std::same_as<std::vector<adjoint_scalar<double>*>>;
    a.for_each_active([](adjoint_scalar<double>&) {});
    a.set_recorded_state(y, t);
  };

// Components [first, first + size) of a System's state.
struct state_range {
  std::size_t first;
  std::size_t size;
};

// A System with a subsystem: a range of its state whose rates read the rest only
// through one input per component, so an explicit step (Cash-Karp or
// Dormand-Prince) can integrate it in substeps of its own under given inputs. A
// step substeps it where the step is past the stepper's stability limit on it,
// where the stepper's error on it set the size of the step before (both unless
// OdeControl::substep_subsystem is off), or where the step spans one of
// `subsystem_knots(from, to)`, the times in (from, to) where its forcing changes
// form; its substeps then land on each knot.
template <typename System>
concept HasSubsystem =
  requires(const System& s, double time,
           const std::vector<typename System::value_type>& y,
           std::vector<typename System::value_type>& out) {
  { s.subsystem() } -> std::same_as<state_range>;
  { s.subsystem_knots(time, time) } -> std::same_as<std::vector<double>>;
  s.subsystem_inputs(out);
  s.subsystem_rates(time, y, y, out);
};

// A System that fills `scales` with one positive, finite multiplier per
// component on its tolerance, for an adaptive step starting at `time`.
template <typename System>
concept ScalesTolerances = requires(const System& s, double time,
                                    std::vector<double>& scales) {
  { s.tolerance_scales(time, scales) } -> std::same_as<void>;
};

// Opt-in domain check. A system may declare
//
//   bool ode_state_valid(const state_type& y) const;
//
// and the adaptive stepper will reject any step landing on a state it refuses,
// shrinking and retrying instead of committing it. Systems that do not declare it
// are unaffected: state_valid() below takes its other branch, so nothing is
// called and nothing costs anything.
//
// The predicate is handed the state *vector* rather than reading the system. The
// stepper's final derivs() does leave the system sitting on y, so either would
// work, but a predicate over a vector is testable without constructing a system
// and is honest about what it is judging.
template <typename System, typename StateType>
concept ChecksState = requires(const System& s, const StateType& y) {
  { s.ode_state_valid(y) } -> std::convertible_to<bool>;
};

template <typename System, typename StateType>
bool state_valid(const System& system, const StateType& y) {
  if constexpr (ChecksState<System, StateType>) {
    return system.ode_state_valid(y);
  } else {
    return true;
  }
}

// Opt-in declaration that the right-hand side does not depend on time. A
// system may declare
//
//   bool ode_autonomous() const;
//
// and the implicit stepper will not form the df/dt term, which it otherwise takes
// by finite difference at the cost of one extra right-hand-side evaluation per
// step. A runtime member rather than a compile-time trait so that one class can
// carry systems of either kind -- the callback system is one class with a flag.
// Systems that do not declare it are treated as time-dependent, which is the
// safe reading.
template <typename System>
class has_autonomous {
  typedef char true_type;
  typedef long false_type;
  template <typename C> static true_type test(decltype(&C::ode_autonomous)) ;
  template <typename C> static false_type test(...);
public:
  enum { value = sizeof(test<System>(0)) == sizeof(true_type) };
};

template <typename System>
typename std::enable_if<has_autonomous<System>::value, bool>::type
is_autonomous(const System& system) {
  return system.ode_autonomous();
}

template <typename System>
typename std::enable_if<!has_autonomous<System>::value, bool>::type
is_autonomous(const System& /* system */) {
  return false;
}

// The recursive interface. Each helper walks a container of elements, threading
// one iterator through them, and is constrained on the one member it calls with
// the iterator it was handed. Constraining the call rather than the element is
// what puts the diagnostic here: an element whose state moves through some other
// scalar's iterator -- a double-typed element reached with an active one -- fails
// the constraint at the call rather than a page of errors inside the loop.
template <typename ForwardIterator>
size_t ode_size(ForwardIterator first, ForwardIterator last) {
  size_t ret = 0;
  while (first != last) {
    ret += first->ode_size();
    ++first;
  }
  return ret;
}

template <typename ForwardIterator>
size_t aux_size(ForwardIterator first, ForwardIterator last) {
  size_t ret = 0;
  while (first != last) {
    ret += first->aux_size();
    ++first;
  }
  return ret;
}

template <typename ForwardIterator, typename It>
  requires requires(std::iter_value_t<ForwardIterator>& e, It it) {
    { e.set_ode_state(it) } -> std::same_as<It>;
  }
It set_ode_state(ForwardIterator first, ForwardIterator last, It it) {
  while (first != last) {
    it = first->set_ode_state(it);
    ++first;
  }
  return it;
}

template <typename ForwardIterator, typename It>
  requires requires(std::iter_value_t<ForwardIterator>& e, It it) {
    { e.ode_state(it) } -> std::same_as<It>;
  }
It ode_state(ForwardIterator first, ForwardIterator last, It it) {
  while (first != last) {
    it = first->ode_state(it);
    ++first;
  }
  return it;
}

template <typename ForwardIterator, typename It>
  requires requires(std::iter_value_t<ForwardIterator>& e, It it) {
    { e.ode_rates(it) } -> std::same_as<It>;
  }
It ode_rates(ForwardIterator first, ForwardIterator last, It it) {
  while (first != last) {
    it = first->ode_rates(it);
    ++first;
  }
  return it;
}

template <typename ForwardIterator, typename It>
  requires requires(std::iter_value_t<ForwardIterator>& e, It it) {
    { e.ode_aux(it) } -> std::same_as<It>;
  }
It ode_aux(ForwardIterator first, ForwardIterator last, It it) {
  while (first != last) {
    it = first->ode_aux(it);
    ++first;
  }
  return it;
}

template <typename ForwardIterator, typename It>
  requires requires(std::iter_value_t<ForwardIterator>& e, It it) {
    { e.set_ode_aux(it) } -> std::same_as<It>;
  }
It set_ode_aux(ForwardIterator first, ForwardIterator last, It it) {
  while (first != last) {
    it = first->set_ode_aux(it);
    ++first;
  }
  return it;
}

template <typename T>
double ode_time(const T& obj) {
  if constexpr (HasOdeTime<T>) {
    return obj.ode_time();
  } else {
    return 0.0;
  }
}

namespace internal {
template <typename T, typename StateType>
void set_ode_state(T& obj, const StateType& y, double time) {
  if constexpr (HasOdeTime<T>) {
    obj.set_ode_state(y.begin(), time);
  } else {
    obj.set_ode_state(y.begin());
  }
}

}

// One rate evaluation: set the state, read the rates. Every stage of a step
// goes through this, and so does every recording a sweep takes.
template <typename T, typename StateType>
void derivs(T& obj, const StateType& y, StateType& dydt,
            const double time) {

  internal::set_ode_state(obj, y, time);
  obj.ode_rates(dydt.begin());
}

// One rate evaluation's use of what it solves for. Opened before the state is
// loaded, because a System's own inner solves can happen inside the load, and
// closed however the evaluation leaves -- which is the whole reason it is a scope
// and not a pair of calls in two places.
//
// Store or load is decided by the CONSTNESS of what it was handed, so there is no
// mode to keep and nothing to hold a stale answer.
//
// `end_solved()` runs from a destructor, so it must not raise: it closes an extent
// and has nothing to fail at.
template <class System, class Values>
struct solved_scope {
  solved_scope(System& system, Values& values, std::size_t stage)
      : system_(system) {
    if constexpr (std::is_const_v<Values>) {
      system_.load_solved(values, stage);
    } else {
      system_.store_solved(values, stage);
    }
  }
  ~solved_scope() { system_.end_solved(); }
  solved_scope(const solved_scope&) = delete;
  solved_scope& operator=(const solved_scope&) = delete;

private:
  System& system_;
};

// One rate evaluation, handed its stage of the step's row to write its inner solves
// into, or to read what an earlier pass wrote; a const row says which.
template <typename T, typename StateType, typename Row>
  requires SolvesForValues<T> &&
           std::same_as<std::remove_const_t<Row>, solved_row<T>>
void derivs(T& obj, const StateType& y, StateType& dydt, const double time,
            Row& row, std::size_t stage) {
  const solved_scope<T, std::remove_reference_t<decltype(row[stage])>> extent{
      obj, row[stage], stage};
  internal::set_ode_state(obj, y, time);
  obj.ode_rates(dydt.begin());
}

// R interface functions - always use std::vector<double>
template <typename T>
std::vector<double> r_derivs(T& obj, const std::vector<double>& y, const double time) {
  std::vector<double> dydt(obj.ode_size());
  derivs(obj, y, dydt, time);
  return dydt;
}

// Two arities rather than one, because R binds each separately: a System with a
// clock is set at a time and one without is not offered one.
template <typename T>
  requires HasOdeTime<T>
void r_set_ode_state(T& obj, const std::vector<double>& y, double time) {
  util::check_length(y.size(), obj.ode_size());
  obj.set_ode_state(y.begin(), time);
}

template <typename T>
  requires (!HasOdeTime<T>)
void r_set_ode_state(T& obj, const std::vector<double>& y) {
  util::check_length(y.size(), obj.ode_size());
  obj.set_ode_state(y.begin());
}

template <typename T>
double r_ode_time(const T& obj) {
  return ode_time(obj);
}

template <typename T>
std::vector<double> r_ode_state(const T& obj) {
  std::vector<double> values(obj.ode_size());
  obj.ode_state(values.begin());
  return values;
}

// Mutable, unlike r_ode_state: a system's ode_rates may compute the rates for
// the state it currently holds rather than return a cached vector.
template <typename T>
std::vector<double> r_ode_rates(T& obj) {
  std::vector<double> dydt(obj.ode_size());
  obj.ode_rates(dydt.begin());
  return dydt;
}

// Put the System on the state the run recorded at `step`. The System reconciles
// itself to that time -- which insertions had happened by then is derived from the
// schedule it was driven by, not handed to it -- and the width it arrives at is
// checked against the width recorded there.
//
// Idempotent, and that is the whole reason it is a load: the System reconciles to
// the step rather than stepping toward it, so arriving twice is arriving once and
// a walk can be run again over the recording it has already walked.
template <class System>
void be_at_step(System& system, std::span<const step_record<System>> rec,
                std::size_t step) {
    if (step >= rec.size()) {
        util::stop("be_at_step: step " +
                   util::to_string(static_cast<int>(step)) +
                   " is outside a recording of " +
                   util::to_string(static_cast<int>(rec.size())) + " steps");
    }
    system.set_recorded_state(rec[step].state, rec[step].time);
    // Named, because a bare length mismatch reads as a caller's error one call
    // away and says nothing about which walk or which step refused.
    if (system.ode_size() != rec[step].state.size()) {
        util::stop("be_at_step: reconciled to " +
                   util::to_string(static_cast<int>(system.ode_size())) +
                   " wide at step " +
                   util::to_string(static_cast<int>(step)) + " against " +
                   util::to_string(static_cast<int>(rec[step].state.size())) +
                   " recorded there");
    }
}

// Apply the insertion the run took at `time`, from the state below it, and read
// back the state it produced: the insertion as a map, so it runs at any scalar
// and a sweep can transpose it.
//
// The System is left holding that wider state -- applying the map is how the state
// is computed -- so a walk rebinds again below this rather than sweeping the width
// below on what this ran on.
//
// A System whose width never changes inserts nothing, so the state passes
// through. That is not a fallback for a System that forgot to declare one -- it
// is what an insertion is for a width that does not move, and it is what lets a
// recording of such a System be walked without it implementing a map it is never
// asked for.
template <class System, class It>
void apply_insertion(System& system, double time, It x,
                     state_type<System>& out) {
  if constexpr (requires { system.apply_insertion(time, x, out); }) {
    system.apply_insertion(time, x, out);
  } else {
    for (std::size_t i = 0; i < out.size(); ++i) {
      out[i] = *x++;
    }
  }
}

template <typename T>
std::vector<double> r_ode_aux(const T& obj) {
  std::vector<double> dydt(obj.aux_size());
  obj.ode_aux(dydt.begin());
  return dydt;
}

template <typename T>
void r_set_ode_aux(T& obj, const std::vector<double>& aux) {
  util::check_length(aux.size(), obj.aux_size());
  obj.set_ode_aux(aux.begin());
}

}
}

#endif
