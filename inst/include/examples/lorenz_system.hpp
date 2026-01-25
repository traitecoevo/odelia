#ifndef LORENZ_SYSTEM_HPP_
#define LORENZ_SYSTEM_HPP_

#include <odelia/ode_solver.hpp>
#include <XAD/XAD.hpp>

using namespace odelia;

// Templated Lorenz system - works with double or XAD active types
template <typename T = double>
class LorenzSystem {
public:
  using value_type = T;  // Expose scalar type for traits
  
  LorenzSystem(T sigma_, T R_, T b_)
    : sigma(sigma_), R(R_), b(b_),
      y0(1.0), y1(1.0), y2(1.0),
      dy0dt(0.0), dy1dt(0.0), dy2dt(0.0),
      time(0.0) {
  }

  // ODE interface
  size_t ode_size() const { return ode_dimension; }

  double ode_time() const { return time; }

  // Iterator methods remain templated on iterator type
  template <typename Iterator>
  Iterator set_ode_state(Iterator it, double time_) {
    time = time_;
    
    y0 = *it++;
    y1 = *it++;
    y2 = *it++;
    
    // Lorenz equations - work with T automatically
    dy0dt = sigma * (y1 - y0);
    dy1dt = R * y0 - y1 - y0 * y2;
    dy2dt = -b * y2 + y0 * y1;
    
    return it;
  }

  template <typename Iterator>
  Iterator ode_state(Iterator it) const {
    *it++ = y0;
    *it++ = y1;
    *it++ = y2;
    return it;
  }

  template <typename Iterator>
  Iterator ode_rates(Iterator it) const {
    *it++ = dy0dt;
    *it++ = dy1dt;
    *it++ = dy2dt;
    return it;
  }

  // Helper to extract double value from any scalar type
  static double to_double(const T& x) {
    // For double: return as-is
    // For XAD: extract value
    if constexpr (std::is_same<T, double>::value) {
      return x;
    } else {
      return xad::value(x);
    }
  }

  std::vector<double> record_step() const {
    std::vector<double> ret;
    ret.reserve(7);
    
    ret.push_back(time);
    ret.push_back(to_double(y0));
    ret.push_back(to_double(y1));
    ret.push_back(to_double(y2));
    ret.push_back(to_double(dy0dt));
    ret.push_back(to_double(dy1dt));
    ret.push_back(to_double(dy2dt));
  
    return ret;
  }

  std::vector<double> pars() const {
    std::vector<double> ret;
    ret.push_back(to_double(sigma));
    ret.push_back(to_double(R));
    ret.push_back(to_double(b));
    return ret;
  }

  void reset() {}

private:
  static const int ode_dimension = 3;

  double time;  // Time always stays double
  T sigma, R, b;  // Parameters can be AD types
  T y0, y1, y2;  // State can be AD types
  T dy0dt, dy1dt, dy2dt;  // Rates can be AD types
};

#endif