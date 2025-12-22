#ifndef LORENZ_SYSTEM_HPP_
#define LORENZ_SYSTEM_HPP_

#include <odelia/ode_solver.hpp>

// define your system using the same namespace as the package 
// makes everything easier to find
namespace odelia {
namespace ode {

class LorenzSystem {
public:
  LorenzSystem(double sigma_, double R_, double b_)
    : sigma(sigma_), R(R_), b(b_),
      y0(1.0), y1(1.0), y2(1.0),
      dy0dt(0.0), dy1dt(0.0), dy2dt(0.0) {
  }

  // ODE interface.
  size_t ode_size() const {return ode_dimension;}

  ode::const_iterator set_ode_state(ode::const_iterator it) {
    y0 = *it++;
    y1 = *it++;
    y2 = *it++;
    
    dy0dt = sigma * (y1 - y0);
    dy1dt = R * y0 - y1 - y0 * y2;
    dy2dt = -b * y2 + y0 * y1;
    
    return it;
  }

  ode::iterator ode_state(ode::iterator it) const {
    *it++ = y0;
    *it++ = y1;
    *it++ = y2;
    return it;
  }

  ode::iterator ode_rates(ode::iterator it) const {
    *it++ = dy0dt;
    *it++ = dy1dt;
    *it++ = dy2dt;
    return it;
  }

  std::vector<double> record_step() const {
    
    std::vector<double> ret;
    
    //todo: add time
    double time = 0 ; // placeholder for time

    ret.push_back(time);
    ret.push_back(y0);
    ret.push_back(y1);
    ret.push_back(y2);
    ret.push_back(dy0dt);
    ret.push_back(dy1dt);
    ret.push_back(dy2dt);
  
    return ret;
  }

  std::vector<double> pars() const {
    std::vector<double> ret;
    ret.push_back(sigma);
    ret.push_back(R);
    ret.push_back(b);
    return ret;
  }

  void reset() {}

  private:
    static const int ode_dimension = 3;

    double sigma, R, b;
    double y0, y1, y2;
    double dy0dt, dy1dt, dy2dt;
  };
}
}
#endif
