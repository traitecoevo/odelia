#ifndef LORENZ_SYSTEM_HPP_
#define LORENZ_SYSTEM_HPP_

#include <odelia/ode_solver.hpp>
#include <XAD/XAD.hpp>
#include <vector>
#include <string>

using namespace odelia;

// Templated Lorenz system - works with double or XAD active types
template <typename T = double>
class LorenzSystem {
public:
  using value_type = T; 
  
  LorenzSystem(T sigma_, T R_, T b_)
    : y0_init(1.0), y1_init(1.0), y2_init(1.0),
      t0(0.0),
      sigma(sigma_), R(R_), b(b_),
      dy0dt(0.0), dy1dt(0.0), dy2dt(0.0) {
    reset();  // initialises state & time
  }

  // rebind names this System on a different scalar, and rebind_from copies its
  // configuration (values only) into that copy. The gradient driver uses them to
  // build the active (double -> AD) version of any System the same way, so a new
  // System gets gradients just by providing these two members. Only values cross,
  // so the copy starts with no tape state; the driver seeds the active inputs after.
  template <class S2> using rebind = LorenzSystem<S2>;

  template <class S2>
  rebind<S2> rebind_from() const {
    // read the parameters and initial state back to plain double (xad::value), build
    // the S2 copy from them, and set its initial state.
    rebind<S2> out(xad::value(sigma), xad::value(R), xad::value(b));
    const double ic[] = {xad::value(y0_init), xad::value(y1_init), xad::value(y2_init)};
    out.set_initial_state(ic, t0);
    return out;
  }

  // ODE interface
  size_t ode_size() const { return ode_dimension; }

  double ode_time() const { return time; }

  double ode_t0() const { return t0; }

  template <typename Iterator>
  Iterator set_ode_state(Iterator it, double time_) {
    time = time_;
    
    y0 = *it++;
    y1 = *it++;
    y2 = *it++;
    
    compute_rates();
    return it;
  }

  void compute_rates() {
    dy0dt = sigma * (y1 - y0);
    dy1dt = R * y0 - y1 - y0 * y2;
    dy2dt = -b * y2 + y0 * y1;
  }

  template <typename Iterator>
  Iterator set_initial_state(Iterator it, double t0_ = 0.0) {
    t0 = t0_;
    y0_init = *it++;
    y1_init = *it++;
    y2_init = *it++;
    return it;
  }

  template <typename Iterator>
  Iterator set_params(Iterator it) {
    sigma = *it++;
    R = *it++;
    b = *it++;
    return it;
  }

  // The differentiable inputs, in the order DifferentiationTargets indexes them:
  // parameters (sigma, R, b) then initial state (y0, y1, y2). The driver seeds a
  // chosen subset active before each solve.
  std::vector<T*> ad_parameters()    { return {&sigma, &R, &b}; }
  std::vector<T*> ad_initial_state() { return {&y0_init, &y1_init, &y2_init}; }

  template <typename Iterator>
  Iterator ode_state(Iterator it) const {
    *it++ = y0;
    *it++ = y1;
    *it++ = y2;
    return it;
  }

  template <typename Iterator>
  Iterator ode_initial_state(Iterator it) const {
    *it++ = y0_init;
    *it++ = y1_init;
    *it++ = y2_init;
    return it;
  }

  template <typename Iterator>
  Iterator ode_rates(Iterator it) const {
    *it++ = dy0dt;
    *it++ = dy1dt;
    *it++ = dy2dt;
    return it;
  }

  std::vector<std::string> record_colnames() const {
    return {"time", "x", "y", "z", "dxdt", "dydt", "dzdt"};
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

  void reset() {
    y0 = y0_init;
    y1 = y1_init;
    y2 = y2_init;
    time = t0;
    compute_rates();
  }

private:
  static const int ode_dimension = 3;

  T y0_init, y1_init, y2_init;
  double t0; 

  T sigma, R, b; 
  T y0, y1, y2; 
  T dy0dt, dy1dt, dy2dt;

  double time;   
};

#endif
