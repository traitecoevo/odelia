// -*-c++-*-
#ifndef ODELIA_ODE_INTERFACE_HPP_
#define ODELIA_ODE_INTERFACE_HPP_

#include <odelia/ode_util.hpp>

#include <concepts>

namespace odelia {
namespace ode {

// Type alias for state vectors based on System's value_type
template<typename System>
using state_type = std::vector<typename System::value_type>;

// Legacy typedefs for R interface (always double)
typedef std::vector<double>::const_iterator const_iterator;
typedef std::vector<double>::iterator       iterator;

// By default, we assume that systems are time homogeneous
template <typename T>
class needs_time {
  typedef char true_type;
  typedef long false_type;
  template <typename C> static true_type test(decltype(&C::ode_time)) ;
  template <typename C> static false_type test(...);
public:
  enum { value = sizeof(test<T>(0)) == sizeof(true_type) };
};

// A System that records the node positions its adaptive solve chose, then replays
// them on a later pass so the active scalar propagates through a schedule that no
// longer moves. The hooks are detected at compile time; an absent hook is a
// zero-cost no-op, so a System that does not record is unaffected. has_recorded_field()
// is required by the concept so a System that provides the hooks but not the query is
// rejected here, at the boundary, rather than deep inside derivs.
//
// Two replay depths. With only positions recorded, the field is recomputed with the
// active scalar on the fixed positions (its derivative flows). With field values
// also recorded, they are reused as fixed doubles (the derivative through the field
// is then zero). has_recorded_field() reports which of the two applies.
template <typename System>
concept Replayable = requires(System s, int stage) {
  s.record_stage(stage);     // per RK stage: record a field value, when values are kept
  s.record_ode_step();       // per accepted ODE step: commit the node positions
  s.replay_step();           // per step on the replay pass: restore this step's record
  { s.has_recorded_field() } -> std::convertible_to<bool>;  // are field values recorded?
};

// The recursive interface
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

template <typename ForwardIterator>
const_iterator set_ode_state(ForwardIterator first, ForwardIterator last,
                             const_iterator it) {
  while (first != last) {
    it = first->set_ode_state(it);
    ++first;
  }
  return it;
}

template <typename ForwardIterator>
iterator ode_state(ForwardIterator first, ForwardIterator last,
                   iterator it) {
  while (first != last) {
    it = first->ode_state(it);
    ++first;
  }
  return it;
}

template <typename ForwardIterator>
iterator ode_rates(ForwardIterator first, ForwardIterator last,
                   iterator it) {
  while (first != last) {
    it = first->ode_rates(it);
    ++first;
  }
  return it;
}

template <typename ForwardIterator>
iterator ode_aux(ForwardIterator first, ForwardIterator last,
                   iterator it) {
  while (first != last) {
    it = first->ode_aux(it);
    ++first;
  }
  return it;
}

template <typename T>
typename std::enable_if<needs_time<T>::value, double>::type
ode_time(const T& obj) {
  return obj.ode_time();
}

template <typename T>
typename std::enable_if<!needs_time<T>::value, double>::type
ode_time(const T& /* obj */) {
  return 0.0;
}

namespace internal {
template <typename T, typename StateType>
typename std::enable_if<needs_time<T>::value, void>::type
set_ode_state(T& obj, const StateType& y, double time) {
  obj.set_ode_state(y.begin(), time);
}

template <typename T, typename StateType>
typename std::enable_if<!needs_time<T>::value, void>::type
set_ode_state(T& obj, const StateType& y, double /* time */) {
  obj.set_ode_state(y.begin());
}

template <typename T, typename StateType>
  requires Replayable<T>
void set_ode_state(T& obj, const StateType& y, int index) {
  obj.set_ode_state(y.begin(), index);
}
}

// primarily for Ode_R - maybe remove
template <typename T, typename StateType>
void derivs(T& obj, const StateType& y, StateType& dydt,
            const double time) {

  internal::set_ode_state(obj, y, time);
  obj.ode_rates(dydt.begin());
}

// ODE stepping. A Replayable System that has recorded field values reads the field
// for this RK stage by index; otherwise it sets state at the current time and
// recomputes (the second branch also covers every non-Replayable System). The choice
// compiles away for a System that isn't Replayable. See the Replayable concept above.
template <typename T, typename StateType>
void derivs(T& obj, const StateType& y, StateType& dydt,
            const double time, const int index) {
  if constexpr (Replayable<T>) {
    if (obj.has_recorded_field()) {
      internal::set_ode_state(obj, y, index);
    } else {
      internal::set_ode_state(obj, y, time);
    }
  } else {
    internal::set_ode_state(obj, y, time);
  }
  obj.ode_rates(dydt.begin());
}

// R interface functions - always use std::vector<double>
template <typename T>
std::vector<double> r_derivs(T& obj, const std::vector<double>& y, const double time) {
  std::vector<double> dydt(obj.ode_size());
  derivs(obj, y, dydt, time);
  return dydt;
}

template <typename T>
typename std::enable_if<needs_time<T>::value, void>::type
r_set_ode_state(T& obj, const std::vector<double>& y, double time) {
  util::check_length(y.size(), obj.ode_size());
  obj.set_ode_state(y.begin(), time);
}

template <typename T>
typename std::enable_if<!needs_time<T>::value, void>::type
r_set_ode_state(T& obj, const std::vector<double>& y) {
  util::check_length(y.size(), obj.ode_size());
  obj.set_ode_state(y.begin());
}

template <typename T>
typename std::enable_if<needs_time<T>::value, double>::type
r_ode_time(const T& obj) {
  return obj.ode_time();
}

template <typename T>
typename std::enable_if<!needs_time<T>::value, double>::type
r_ode_time(const T& /* obj */) {
  return 0.0;
}

template <typename T>
std::vector<double> r_ode_state(const T& obj) {
  std::vector<double> values(obj.ode_size());
  obj.ode_state(values.begin());
  return values;
}

template <typename T>
std::vector<double> r_ode_rates(const T& obj) {
  std::vector<double> dydt(obj.ode_size());
  obj.ode_rates(dydt.begin());
  return dydt;
}

template <typename T>
std::vector<double> r_ode_aux(const T& obj) {
  std::vector<double> dydt(obj.aux_size());
  obj.ode_aux(dydt.begin());
  return dydt;
}

}
}

#endif
