#ifndef ODELIA_ODE_CONTROL_HPP_
#define ODELIA_ODE_CONTROL_HPP_

#include <vector>
#include <cstddef>
#include <odelia/ode_util.hpp>

namespace odelia {
namespace ode {

struct OdeControl {
  typedef std::vector<double> state_type;

  OdeControl() : OdeControl(1e-8, 1e-8, 1.0, 0.0,
                            1e-8, 10.0, 1e-6)
  {
  }
  
  OdeControl(double tol_abs_, double tol_rel_,
             double a_y_, double a_dydt_,
             double step_size_min_, double step_size_max_,
             double step_size_initial_)
  {
    set_controls(tol_abs_, tol_rel_, a_y_, a_dydt_,
                 step_size_min_, step_size_max_, step_size_initial_);
  }

  void set_controls(double tol_abs_, double tol_rel_,
                    double a_y_, double a_dydt_,
                    double step_size_min_, double step_size_max_,
                    double step_size_initial_) {
    set_tol_rel(tol_rel_);
    set_tol_abs(tol_abs_);
    set_a_y(a_y_);
    set_a_dydt(a_dydt_);
    set_step_size_min(step_size_min_);
    set_step_size_max(step_size_max_);
    set_step_size_initial(step_size_initial_);
    last_step_size_shrank = false;
  }

  std::vector<double> get_controls() const {
    return {tol_abs, tol_rel, a_y, a_dydt, 
            step_size_min, step_size_max, step_size_initial};
  }

  void set_tol_abs(double value) { tol_abs = value; }
  double get_tol_abs() const { return tol_abs; }

  void set_tol_rel(double value) { tol_rel = value; }
  double get_tol_rel() const { return tol_rel; }

  void set_a_y(double value) { a_y = value; }
  double get_a_y() const { return a_y; }

  void set_a_dydt(double value) { a_dydt = value; }
  double get_a_dydt() const { return a_dydt; }

  void set_step_size_min(double value) { step_size_min = value; }
  double get_step_size_min() const { return step_size_min; }

  void set_step_size_max(double value) { step_size_max = value; }
  double get_step_size_max() const { return step_size_max; }

  void set_step_size_initial(double value) { step_size_initial = value; }
  double get_step_size_initial() const { return step_size_initial; }

  double adjust_step_size(size_t dim, size_t ord,
                          double step_size,
                          const state_type &y,
                          const state_type &yerr,
                          const state_type &dydt)
  {
    double rmax = std::numeric_limits<double>::min();
    const double S = 0.9;

    for (size_t i = 0; i < dim; i++)
    {
      const double D0 = errlevel(y[i], dydt[i], step_size);
      using std::abs;
      const double r = abs(yerr[i]) / abs(D0);
      rmax = std::max(r, rmax);
    }

    if (rmax > 1.1)
    {
      // decrease step, no more than factor of 5
      double r = S / pow(rmax, 1.0 / ord);
      if (r < 0.2)
      {
        r = 0.2;
      }
      double new_step = step_size * r;
      if (new_step < step_size_min)
      {
        new_step = step_size_min;
      }

      if (new_step < step_size)
      {
        step_size = new_step;
        last_step_size_shrank = true;
      }
    }
    else if (rmax < 0.5)
    {
      // increase step, no more than factor of 5
      double r = S / pow(rmax, 1.0 / (ord + 1.0));
      if (r > 5.0)
      {
        r = 5.0;
      }
      else if (r < 1.0)
      {
        r = 1.0;
      }
      step_size *= r;
      if (step_size > step_size_max)
      {
        step_size = step_size_max;
      }
      last_step_size_shrank = false;
    }
    else
    {
      last_step_size_shrank = false;
    }

    return step_size;
  }

  double errlevel(double y, double dydt, double h) const
  {
    using std::abs;
    const double errlev = tol_rel * (a_y * abs(y) +
                                    a_dydt * abs(h * dydt)) +
                          tol_abs;
    if (errlev <= 0.0)
    {
      util::stop("errlev <= zero");
    }
    return errlev;
  }

  bool step_size_shrank() const
  {
    return last_step_size_shrank;
  }
  
  double tol_abs, tol_rel, a_y, a_dydt;
  double step_size_min, step_size_max, step_size_initial;
  bool last_step_size_shrank;
};

}
}
#endif