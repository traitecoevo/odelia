#ifndef LEAF_THERMAL_SYSTEM_HPP_
#define LEAF_THERMAL_SYSTEM_HPP_

#include <odelia/ode_solver.hpp>
#include <odelia/drivers.hpp>
#include <XAD/XAD.hpp>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace odelia;

// Parameter struct - always doubles, used to initialise the system
struct LeafThermalPars {
  double k_H = 0.5;       // 1/h
  double g_tr_max = 1.0;  // °C/h cooling
  double m_tr = 0.5;      // steepness
  double T_tr_mid = 30.0; // °C
};

// Templated system - works with double or XAD active types
template <typename T = double>
class LeafThermalSystem {
public:
  using value_type = T;

  LeafThermalSystem(const LeafThermalPars &pars_, const drivers::Drivers &drivers_)
      : pars(pars_),
        k_H(pars_.k_H),
        g_tr_max(pars_.g_tr_max),
        m_tr(pars_.m_tr),
        T_tr_mid(pars_.T_tr_mid),
        T_LC_init(0.0),
        t0(0.0),
        dT_LC(0.0),
        T_air(0.0),
        S_tr(0.0) {
    initialize_drivers(drivers_);
    reset();
  }

  void initialize_drivers(const drivers::Drivers &drv) {
    drivers = drv;
    temperature_fn = drv.get_function_ptr("temperature");
    if (!temperature_fn)
      throw std::runtime_error("Missing driver 'temperature' for LeafThermalSystem");
  }

  // ODE interface
  size_t ode_size() const { return ode_dimension; }
  double ode_time() const { return time; }
  double ode_t0() const { return t0; }

  template <typename Iterator>
  Iterator set_ode_state(Iterator it, double time_) {
    T_LC = *it++;
    time = time_;
    set_drivers();
    compute_ode_rates();
    return it;
  }

  // Set the initial state (reset point) - no tape registration
  template <typename Iterator>
  Iterator set_initial_state(Iterator it, double t0_ = 0.0) {
    t0 = t0_;
    T_LC_init = *it++;
    return it;
  }

  // Set the initial state and register on tape for AD gradient computation
  template <typename Tape, typename Iterator>
  std::vector<T*> set_initial_state(Tape& tape, Iterator it, double t0_) {
    t0 = t0_;
    T_LC_init = *it++;
    tape.registerInput(T_LC_init);
    return {&T_LC_init};
  }

  // Set parameters - no tape registration
  template <typename Iterator>
  Iterator set_params(Iterator it) {
    k_H = *it++;
    g_tr_max = *it++;
    m_tr = *it++;
    T_tr_mid = *it++;
    return it;
  }

  // Set parameters and register on tape for AD gradient computation
  template <typename Tape, typename Iterator>
  std::vector<T*> set_params(Tape& tape, Iterator it) {
    k_H = *it++;
    g_tr_max = *it++;
    m_tr = *it++;
    T_tr_mid = *it++;
    tape.registerInput(k_H);
    tape.registerInput(g_tr_max);
    tape.registerInput(m_tr);
    tape.registerInput(T_tr_mid);
    return {&k_H, &g_tr_max, &m_tr, &T_tr_mid};
  }

void set_drivers() {
    T_air = temperature_fn->evaluate(time);
  }

  T logistic_raw(T z, T slope) {
    return 1.0 / (1.0 + exp(slope * z));
  }

  void compute_ode_rates() {
    S_tr = logistic_raw(T_LC - T_tr_mid, -m_tr);
    dT_LC = k_H * (T_air - T_LC) - g_tr_max * S_tr;
  }

  template <typename Iterator>
  Iterator ode_state(Iterator it) const {
    *it++ = T_LC;
    return it;
  }

  template <typename Iterator>
  Iterator ode_initial_state(Iterator it) const {
    *it++ = T_LC_init;
    return it;
  }

  template <typename Iterator>
  Iterator ode_rates(Iterator it) const {
    *it++ = dT_LC;
    return it;
  }

  std::vector<double> get_current_drivers() const {
    std::vector<double> ret;
    ret.push_back(temperature_fn->evaluate(time));
    return ret;
  }

  std::vector<std::string> record_colnames() const {
    return {"time", "T_LC", "T_air", "dT_LC", "S_tr"};
  }

  std::vector<double> record_step() const {
    std::vector<double> ret;
    ret.push_back(time);
    ret.push_back(xad::value(T_LC));
    ret.push_back(T_air);
    ret.push_back(xad::value(dT_LC));
    ret.push_back(xad::value(S_tr));
    return ret;
  }

  std::vector<double> get_pars() const {
    std::vector<double> ret;
    ret.push_back(xad::value(k_H));
    ret.push_back(xad::value(g_tr_max));
    ret.push_back(xad::value(m_tr));
    ret.push_back(xad::value(T_tr_mid));
    return ret;
  }

  void reset() {
    T_LC = T_LC_init;
    time = t0;
    set_drivers();
    compute_ode_rates();
  }

private:
  static const int ode_dimension = 1;

  LeafThermalPars pars;

  // Parameters (can be AD types)
  T k_H, g_tr_max, m_tr, T_tr_mid;

  // Initial state and time (reset point)
  T T_LC_init;
  double t0;

  // Current state and rates (can be AD types)
  T T_LC;
  T dT_LC;

  // Time and drivers (always double)
  double time;
  drivers::Drivers drivers;
  const drivers::Function *temperature_fn;

  // Auxiliary
  double T_air;
  T S_tr;
};

#endif