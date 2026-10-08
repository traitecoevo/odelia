#ifndef ODELIA_DRIVERS_H
#define ODELIA_DRIVERS_H

#include <string>
#include <unordered_map>
#include <vector>
#include <odelia/interpolator.hpp>

namespace odelia {
namespace drivers {

// How a driver given as values alone gets its knot slopes. `natural` is the
// natural cubic spline, the default and what odelia has always read. `monotone`
// (Fritsch-Carlson) keeps every span inside the two values bracketing it, so a
// series that is never negative -- rainfall -- is never read negative between its
// points; a natural spline beside a wet day can dip below zero.
enum class Slopes { natural, monotone };

class Function
{
public:
  Function() = default;

  Function(std::vector<double> const &x, std::vector<double> const &y,
           Slopes slopes = Slopes::natural)
  {
    if (slopes == Slopes::monotone)
    {
      variable.init(x, y, interpolator::monotone_slopes<double>(x, y));
    }
    else
    {
      variable.init(x, y);
    }
    variable.set_extrapolate(false);
    is_variable = true;
  }

  Function(double k)
  {
    constant = k;
    is_variable = false;
  }

  double evaluate(double u) const
  {
    if (is_variable)
    {
      return variable.eval(u);
    }
    else
    {
      return constant;
    }
  }

  std::vector<double> evaluate_range(const std::vector<double> &u) const
  {
    if (is_variable)
    {
      return variable.r_eval(u);
    }
    else
    {
      return std::vector<double>(u.size(), constant);
    }
  }

  void set_extrapolate(bool extrapolate)
  {
    variable.set_extrapolate(extrapolate);
  }

private:
  interpolator::Interpolator variable;
  double constant;
  bool is_variable;
};

class Drivers {
  
public:
  // this will override any previously defined drivers with the same name
  void set_constant(std::string driver_name, double k) {
   
    if (drivers.find(driver_name) != drivers.end()) 
    {
      drivers.erase(driver_name);
    }
    drivers.insert({driver_name, drivers::Function(k)});
  }

  // initialise spline of driver with x, y control points; `slopes` says how the
  // knot slopes are chosen (see Slopes)
  void set_variable(std::string driver_name, std::vector<double> const &x,
                    std::vector<double> const &y,
                    Slopes slopes = Slopes::natural) {
    if (drivers.find(driver_name) != drivers.end())
    {
      drivers.erase(driver_name);
    }
    drivers.insert({driver_name, drivers::Function(x, y, slopes)});
  }

  void set_extrapolate(std::string driver_name, bool extrapolate) {
    drivers.at(driver_name).set_extrapolate(extrapolate);
  }

  // evaluate/query interpolated spline for driver at point u, return s(x), where s is interpolated function
  double evaluate(std::string driver_name, double x) const {

    return drivers.at(driver_name).evaluate(x);
  }

  // evaluate/query interpolated spline for driver at vector of points, return vector of values
  std::vector<double> evaluate_range(std::string driver_name, std::vector<double> x) const {
    return drivers.at(driver_name).evaluate_range(x);
  }

  // returns the name of each active driver - useful for R output
  std::vector<std::string> get_names() const {
    auto ret = std::vector<std::string>();
    for (auto const &driver: drivers) {
      ret.push_back(driver.first);
    }
    return ret;
  }

  void clear() {
    drivers.clear();
  }

  // Get pointer to Function for repeated use (returns nullptr if not found)
  const Function *get_function_ptr(const std::string &driver_name) const
  {
    auto it = drivers.find(driver_name);
    if (it != drivers.end())
    {
      return &(it->second);
    }
    return nullptr;
  }

private:
  std::unordered_map <std::string, drivers::Function> drivers;
};
}

}

#endif //ODELIA_DRIVERS_H
