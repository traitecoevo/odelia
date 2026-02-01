// -*-c++-*-
#ifndef ODELIA_ODE_INTERFACE_HPP_
#define ODELIA_ODE_INTERFACE_HPP_

#include <odelia/ode_util.hpp>

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
      internal::set_ode_state(obj, y, index);
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
