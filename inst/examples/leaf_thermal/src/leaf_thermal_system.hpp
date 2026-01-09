#include <odelia/ode_solver.hpp>
#include <odelia/drivers.hpp>

#include "leaf_thermal_pars.hpp"
#include <vector>
#include <cmath>
#include <algorithm>

namespace odelia
{
namespace ode
{

class LeafThermalSystem {
  public:
    // Constructor: pass parameter struct
    LeafThermalSystem(const LeafThermalPars &pars_)
        : pars(pars_),
          T_LC(0.0),
          dT_LC(0.0),
          time(0.0),
          // auxilliary variables
          T_air(0.0),
          S_tr(0.0)
    {
    }

    size_t ode_size() const { return ode_dimension; }

    double ode_time() const { return time;}

    ode::const_iterator set_ode_state(ode::const_iterator it, double time_)
    {
      // Unpack state
      T_LC = *it++;
      time = time_;

      // Compute rates
      compute_ode_rates();

      return it;
    }
    
    double logistic_raw(double z, double slope)
    {
      return 1.0 / (1.0 + std::exp(slope * z));
    }

    void compute_ode_rates()
    {
      const double pi = 3.14159265358979323846;

      // Air temperature sinusoid
      double Tmean = 32.0; // °C
      double Tamp = 6.0;   // °C amplitude
      double tpeak = 15.0; // hour of daily peak

      // 1. Air temperature forcing
      T_air = Tmean + Tamp *
                                std::sin(2.0 * pi * (time - tpeak) / 24.0);
      // 2. Transpiration rate S_tr (logistic function)
      S_tr = logistic_raw(T_LC - pars.T_tr_mid, -pars.m_tr);
      dT_LC = pars.k_H * (T_air - T_LC) - pars.g_tr_max * S_tr;
    }

    ode::iterator ode_state(ode::iterator it) const
    {
      *it++ = T_LC;

      return it;
    }

    ode::iterator ode_rates(ode::iterator it) const
    {
      *it++ = dT_LC;

      return it;
    }

    std::vector<std::string> record_colnames() const {
      return {
          "time",
          "T_LC",
          "T_air",
          "dT_LC",
      };
    }

  std::vector<double> record_step() const
    {

      std::vector<double> ret;

      ret.push_back(time);
      ret.push_back(T_LC);
      ret.push_back(T_air);
      ret.push_back(dT_LC);

      return ret;
    }

    // Optional: return parameters as a vector
    std::vector<double> get_pars() const
    {
      const double *p = reinterpret_cast<const double *>(&pars);
      return std::vector<double>(p, p + sizeof(LeafThermalPars) / sizeof(double));
    }

    void reset() {}

    // Diagnostics if needed
    double get_T_air() const { return T_air; }

    LeafThermalPars pars;

    drivers::Drivers drivers;

  private:
    static const int ode_dimension = 1;

    // State
    double T_LC;

    // Rates
    double dT_LC;

    // Time for forcing
    double time;

    // Auxiliary variables
    double T_air;
    double S_tr;
    };
} // namespace ode
} // namespace odelia
