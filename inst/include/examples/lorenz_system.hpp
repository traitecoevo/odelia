#ifndef LORENZ_SYSTEM_HPP_
#define LORENZ_SYSTEM_HPP_

#include <odelia/ode_solver.hpp>
#include <XAD/XAD.hpp>

using namespace odelia;

// Templated Lorenz system - works with double or XAD active types
template <typename T = double>
class LorenzSystem {
public:
  using value_type = T; 
  
  LorenzSystem(T sigma_, T R_, T b_)
    : sigma(sigma_), R(R_), b(b_),
      y0(1.0), y1(1.0), y2(1.0),
      dy0dt(0.0), dy1dt(0.0), dy2dt(0.0),
      time(0.0) {
  }

  // ODE interface
  size_t ode_size() const { return ode_dimension; }

  double ode_time() const { return time; }

  template <typename Iterator>
  Iterator set_ode_state(Iterator it, double time_) {
    time = time_;
    
    y0 = *it++;
    y1 = *it++;
    y2 = *it++;
    
    compute_rates();
    return it;
  }

  // Registers inputs, returns pointers for AD gradient computation
  template <typename Tape, typename Iterator>
  std::vector<T*> set_ode_state(Tape& tape, Iterator it, double time_) {
    time = time_;
    
    y0 = *it++;
    y1 = *it++;
    y2 = *it++;
    
    compute_rates();
    
    tape.registerInput(y0);
    tape.registerInput(y1);
    tape.registerInput(y2);
    return {&y0, &y1, &y2};
  }

  void compute_rates() {
    dy0dt = sigma * (y1 - y0);
    dy1dt = R * y0 - y1 - y0 * y2;
    dy2dt = -b * y2 + y0 * y1;
  }

  template <typename Iterator>
  Iterator set_params(Iterator it) {
    sigma = *it++;
    R = *it++;
    b = *it++;
    return it;
  }

  // Registers inputs, returns pointers for AD gradient computation
  template <typename Tape, typename Iterator>
  std::vector<T*> set_params(Tape& tape, Iterator it) {
    sigma = *it++;
    R = *it++;
    b = *it++;
    tape.registerInput(sigma);
    tape.registerInput(R);
    tape.registerInput(b);
    return {&sigma, &R, &b};
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

  std::vector<double> record_step() const {
    std::vector<double> ret;
    ret.reserve(7);
    
    ret.push_back(time);
    ret.push_back(xad::value(y0));
    ret.push_back(xad::value(y1));
    ret.push_back(xad::value(y2));
    ret.push_back(xad::value(dy0dt));
    ret.push_back(xad::value(dy1dt));
    ret.push_back(xad::value(dy2dt));
    
    return ret;
  }

  std::vector<double> pars() const {
    std::vector<double> ret;
    ret.push_back(xad::value(sigma));
    ret.push_back(xad::value(R));
    ret.push_back(xad::value(b));
    return ret;
  }

  void reset() {}

private:
  static const int ode_dimension = 3;

  double time; 
  T sigma, R, b; 
  T y0, y1, y2; 
  T dy0dt, dy1dt, dy2dt; 
};

#endif