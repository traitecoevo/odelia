#ifndef LEAF_THERMAL_SYSTEM_HPP_
#define LEAF_THERMAL_SYSTEM_HPP_

#include <odelia/ode_solver.hpp>
#include <odelia/drivers.hpp>
#include <XAD/XAD.hpp>
#include <memory>
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

  // rebind + rebind_from copy this System's configuration onto another scalar, so
  // the gradient driver can build its active (AD) version generically (see
  // LorenzSystem for the pattern). Only values cross; the driver seeds the active
  // inputs afterwards.

  template <class S2>
  LeafThermalSystem<S2> rebind_from() const {
    const LeafThermalPars p{xad::value(k_H), xad::value(g_tr_max),
                            xad::value(m_tr), xad::value(T_tr_mid)};
    LeafThermalSystem<S2> out(p, *drivers);
    const double ic = xad::value(T_LC_init);
    out.set_initial_state(&ic, t0);
    return out;
  }

  void initialize_drivers(const drivers::Drivers &drv) {
    drivers = std::make_shared<const drivers::Drivers>(drv);
    temperature_fn = drivers->get_function_ptr("temperature");
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

  // Stand on a state a run recorded, for a reverse sweep: the drivers are read
  // at the recorded time. The width never changes, so this is set_ode_state.
  void set_recorded_state(const std::vector<T>& y, double time_) {
    set_ode_state(y.begin(), time_);
  }

  // Set the initial state (reset point) - no tape registration
  template <typename Iterator>
  Iterator set_initial_state(Iterator it, double t0_ = 0.0) {
    t0 = t0_;
    T_LC_init = *it++;
    return it;
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

  // The parameters a pass can seed active, in the order it indexes them.
  std::vector<T*> ad_parameters() { return {&k_H, &g_tr_max, &m_tr, &T_tr_mid}; }

  // Every member carrying the scalar, for a reverse sweep. ⚠️ A member left out
  // here silently contributes nothing to the gradient; T_air is a driver value
  // and stays double. Both overloads, because visit_active passes over a const
  // object whose walk is non-const without saying so.
  template <class F>
  void for_each_active(F&& f) {
    f(k_H); f(g_tr_max); f(m_tr); f(T_tr_mid);
    f(T_LC_init); f(T_LC); f(dT_LC); f(S_tr);
  }
  template <class F>
  void for_each_active(F&& f) const {
    f(k_H); f(g_tr_max); f(m_tr); f(T_tr_mid);
    f(T_LC_init); f(T_LC); f(dT_LC); f(S_tr);
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
  // ⚠️ SHARED, not copied, because temperature_fn points into it. A System copied
  // with its own Drivers kept a pointer into the source's, which die with the
  // source -- and R frees the system a solver was built from as soon as nothing
  // holds it, after which every read was of freed memory. Shared, every copy's
  // pointer is into one Drivers that lives as long as any copy does. Re-deriving
  // the pointer after each copy instead measured 2x slower: the solver copies the
  // System often.
  std::shared_ptr<const drivers::Drivers> drivers;
  const drivers::Function *temperature_fn = nullptr;

  // Auxiliary
  double T_air;
  T S_tr;
};

#endif