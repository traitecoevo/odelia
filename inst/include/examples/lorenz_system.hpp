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
  
  // Every scalar's LorenzSystem is one class, so assign_from reaches the source's
  // members.
  template <typename> friend class LorenzSystem;

  LorenzSystem(T sigma_ = T(0.0), T R_ = T(0.0), T b_ = T(0.0))
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

  // The one map: the parameters and the initial state, read back to plain double
  // (xad::value) so only values cross. rebind_from is a line over it.
  template <class S1>
  void assign_from(const LorenzSystem<S1>& src) {
    sigma = T(xad::value(src.sigma));
    R     = T(xad::value(src.R));
    b     = T(xad::value(src.b));
    const double ic[] = {xad::value(src.y0_init), xad::value(src.y1_init),
                         xad::value(src.y2_init)};
    set_initial_state(ic, src.t0);
  }

  template <class S2>
  LorenzSystem<S2> rebind_from() const {
    LorenzSystem<S2> out;
    out.assign_from(*this);
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

  // Stand on a state a run recorded, for a reverse sweep. The width never
  // changes, so this is set_ode_state.
  void set_recorded_state(const std::vector<T>& y, double time_) {
    set_ode_state(y.begin(), time_);
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

  // The parameters a pass can seed active, in the order it indexes them.
  std::vector<T*> ad_parameters() { return {&sigma, &R, &b}; }

  // Everything here that carries the scalar, which is what a walk holding this
  // System across recordings hands back before it clears the tape.
  template <class F>
  void for_each_active(F&& f) {
    f(sigma); f(R); f(b);
    f(y0); f(y1); f(y2);
    f(dy0dt); f(dy1dt); f(dy2dt);
    f(y0_init); f(y1_init); f(y2_init);
  }

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
