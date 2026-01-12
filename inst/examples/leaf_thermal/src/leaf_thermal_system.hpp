#include <odelia/ode_solver.hpp>
#include <odelia/drivers.hpp>

#include "leaf_thermal_pars.hpp"
#include <vector>
#include <cmath>
#include <algorithm>

using namespace odelia;

class LeafThermalSystem {
  public:
    // Constructor: pass parameter struct and drivers
    LeafThermalSystem(const LeafThermalPars &pars_, const drivers::Drivers &drivers_)
        : pars(pars_),
          // Initialized to zero
          // These will be reset in set_ode_state
          T_LC(0.0),
          dT_LC(0.0),
          time(0.0),
          T_air(0.0),
          S_tr(0.0) {
      initialize_drivers(drivers_);
    }

    // Initialize drivers 
    // - Drivers must have been set up externally and passed in here
    // - The data is an unordered_map of driver name to driver Function
    // - So we can look up by name
    // - We can cache pointers to specific functions for speed
    // - This saves repeated lookups during simulation
    void initialize_drivers(const drivers::Drivers &drv) {
      drivers = drv;

      // Cache at setup time
      temperature_fn = drv.get_function_ptr("temperature");
      if (!temperature_fn)
        throw std::runtime_error("Missing driver 'temperature' for LeafThermalSystem");
    }

    size_t ode_size() const { return ode_dimension; }

    double ode_time() const { return time;}

    ode::const_iterator set_ode_state(ode::const_iterator it, double time_) {
      // Unpack state
      T_LC = *it++;
      time = time_;

      // Set drivers based on current time
      set_drivers();
 
      // Compute rates at this time and state
      compute_ode_rates();

      return it;
    }

    // Set drivers at current time
    void set_drivers() {
      T_air = temperature_fn->evaluate(time);
    }

    double logistic_raw(double z, double slope) {
      return 1.0 / (1.0 + std::exp(slope * z));
    }

    // Compute ODE rates based on current state and drivers
    // - Called from set_ode_state after drivers are updated
    // - This is where the model equations are implemented
    void compute_ode_rates() {

      // Transpiration rate varies with temperature 
      // - logistic function
      S_tr = logistic_raw(T_LC - pars.T_tr_mid, -pars.m_tr);
      // Change in leaf temperature
      dT_LC = pars.k_H * (T_air - T_LC) - pars.g_tr_max * S_tr;
    }

    ode::iterator ode_state(ode::iterator it) const {
      *it++ = T_LC;

      return it;
    }

    ode::iterator ode_rates(ode::iterator it) const {
      *it++ = dT_LC;

      return it;
    }

    std::vector<double> get_current_drivers() const {
      std::vector<double> ret;
      ret.push_back(temperature_fn->evaluate(time));
      return ret;
    }

    std::vector<std::string> record_colnames() const {
      return {
          "time",
          "T_LC",
          "T_air",
          "dT_LC",
          "S_tr"
      };
    }

  // Optional: compile a vector of variables to record at each step
  std::vector<double> record_step() const {

      std::vector<double> ret;

      ret.push_back(time);
      ret.push_back(T_LC);
      ret.push_back(T_air);
      ret.push_back(dT_LC);
      ret.push_back(S_tr);

      return ret;
    }

    // Optional: return parameters as a vector
    std::vector<double> get_pars() const {
      const double *p = reinterpret_cast<const double *>(&pars);
      return std::vector<double>(p, p + sizeof(LeafThermalPars) / sizeof(double));
    }

    void reset() {}

  private:
    static const int ode_dimension = 1;
    
    // Parameters
    LeafThermalPars pars;

    // State
    double T_LC;

    // Rates
    double dT_LC;

    // Time and drivers for forcing
    double time;
    // Drivers object - unordered map of driver name to Function
    drivers::Drivers drivers;
    // Cached pointer to specific driver (for speed)
    const drivers::Function *temperature_fn; 

    // Auxiliary variables
    double T_air;
    double S_tr;
};
