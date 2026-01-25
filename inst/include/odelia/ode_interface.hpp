// -*-c++-*-
#ifndef ODELIA_ODE_INTERFACE_HPP_
#define ODELIA_ODE_INTERFACE_HPP_

#include <odelia/ode_util.hpp>

namespace odelia {
namespace ode {

// 
// System Traits: Extract scalar type from ODE systems
// 
// 
// Primary template - default for non-templated systems (backward compatibility)
template<typename System, typename = void>
struct system_traits {
    using value_type = double;
    using state_type = std::vector<value_type>;
};

// Specialization for systems that define value_type (templated systems)
// Uses SFINAE to detect if System::value_type exists
template<typename System>
struct system_traits<System, 
    typename std::enable_if<
        !std::is_same<typename System::value_type, void>::value
    >::type> {
    using value_type = typename System::value_type;
    using state_type = std::vector<value_type>;
};

// 
// Legacy typedefs (kept for backward compatibility with existing code)
//
// These are utilities designed to make it more pleasant to work with
// ode objects.  The first in-place functions work with the
// boost::odeint interface, and the types are just to reduce typing.
typedef std::vector<double>        state_type;
typedef state_type::const_iterator const_iterator;
typedef state_type::iterator       iterator;

// By default, we assume that systems are time homogeneous; systems
// that provide an `ode_time` function will be treated differently.
template <typename T>
class needs_time {
  typedef char true_type;
  typedef long false_type;
  template <typename C> static true_type test(decltype(&C::ode_time)) ;
  template <typename C> static false_type test(...);
public:
  enum { value = sizeof(test<T>(0)) == sizeof(true_type) };
};

// Have a special case where we want to store and reuse Patch state at each RK45 step
template <typename System>
class has_cache {
  typedef char true_type;
  typedef long false_type;
  template <typename C> static true_type test(decltype(&C::cache_RK45_step)) ;
  template <typename C> static false_type test(...);
public:
  enum { value = sizeof(test<System>(0)) == sizeof(true_type) };
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
typename std::enable_if<has_cache<T>::value, void>::type
set_ode_state(T& obj, const StateType& y, int index) {
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

// for ODE stepping
template <typename T, typename StateType>
typename std::enable_if<!has_cache<T>::value, void>::type
derivs(T& obj, const StateType& y, StateType& dydt,
            const double time, const int /* index */) {

  internal::set_ode_state(obj, y, time);
  obj.ode_rates(dydt.begin());
}

// for mutants or ODE stepping
template <typename T, typename StateType>
typename std::enable_if<has_cache<T>::value, void>::type
derivs(T& obj, const StateType& y, StateType& dydt,
            const double time, const int index) {

    if(obj.use_cached_environment) {
      internal::set_ode_state(obj, y, index); // only works for patches
    } else {
      internal::set_ode_state(obj, y, time);
    }
  
  obj.ode_rates(dydt.begin());
}

template <typename T>
state_type r_derivs(T& obj, const state_type& y, const double time) {
  state_type dydt(obj.ode_size());
  derivs(obj, y, dydt, time);
  return dydt;
}

// These out-of-place versions are useful for interfacing with R.
template <typename T>
typename std::enable_if<needs_time<T>::value, void>::type
r_set_ode_state(T& obj, const state_type& y, double time) {
  util::check_length(y.size(), obj.ode_size());
  obj.set_ode_state(y.begin(), time);
}

template <typename T>
typename std::enable_if<!needs_time<T>::value, void>::type
r_set_ode_state(T& obj, const state_type& y) {
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
state_type r_ode_state(const T& obj) {
  state_type values(obj.ode_size());
  obj.ode_state(values.begin());
  return values;
}

template <typename T>
state_type r_ode_rates(const T& obj) {
  state_type dydt(obj.ode_size());
  obj.ode_rates(dydt.begin());
  return dydt;
}

template <typename T>
state_type r_ode_aux(const T& obj) {
  state_type dydt(obj.aux_size());
  obj.ode_aux(dydt.begin());
  return dydt;
}

}
}

#endif