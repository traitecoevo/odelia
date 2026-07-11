/* Generic ODE interface for odelia package
 * 
 * This file provides generic functions that can be used with any ODE system.
 * It includes:
 * - Drivers: External forcing data (constant or time-varying)
 * - OdeControl: Solver settings and tolerances
 * - Generic Solver functions: Templated implementations that work with any System type
 * 
 * These functions are system-agnostic and can be used by any model
 * (Lorenz, Leaf Thermal, ATLAS, etc.)
 * 
 * System-specific interfaces (lorenz_interface.cpp, leaf_thermal_interface.cpp)
 * call these templated functions with their specific System types.
 */

#include <Rcpp.h>
#include <XAD/XAD.hpp>
#include <odelia/drivers.hpp>
#include <odelia/ode_control.hpp>
#include <odelia/ode_solver.hpp>
#include <odelia/gradient.hpp>
#include <odelia/rcpp_interface_helpers.hpp>

using namespace Rcpp;
using namespace odelia;

namespace odelia {
namespace solver {
template<typename T>
inline Rcpp::XPtr<ode::Solver<T>> get_solver(SEXP xp);
}  // namespace solver
}  // namespace odelia

//-------------------------------------------------------------------------
// OdeControl interface

// Helper to get OdeControl pointer
inline Rcpp::XPtr<ode::OdeControl> get_OdeControl(SEXP xp) {
  return Rcpp::XPtr<ode::OdeControl>(xp);
}

// [[Rcpp::export]]
SEXP OdeControl_new() {
  return Rcpp::XPtr<ode::OdeControl>(new ode::OdeControl(), true);
}

// [[Rcpp::export]]
Rcpp::List OdeControl_get_controls(SEXP control_xp) {
  auto ctrl = get_OdeControl(control_xp);
  return Rcpp::List::create(
    Rcpp::Named("tol_abs") = ctrl->tol_abs,
    Rcpp::Named("tol_rel") = ctrl->tol_rel,
    Rcpp::Named("a_y") = ctrl->a_y,
    Rcpp::Named("a_dydt") = ctrl->a_dydt,
    Rcpp::Named("step_size_min") = ctrl->step_size_min,
    Rcpp::Named("step_size_max") = ctrl->step_size_max,
    Rcpp::Named("step_size_initial") = ctrl->step_size_initial
  );
}

// [[Rcpp::export]]
void OdeControl_set_controls(SEXP control_xp,
                             double tol_abs, double tol_rel,
                             double a_y, double a_dydt,
                             double step_size_min, double step_size_max,
                             double step_size_initial) {
  auto ctrl = get_OdeControl(control_xp);
  ctrl->tol_abs = tol_abs;
  ctrl->tol_rel = tol_rel;
  ctrl->a_y = a_y;
  ctrl->a_dydt = a_dydt;
  ctrl->step_size_min = step_size_min;
  ctrl->step_size_max = step_size_max;
  ctrl->step_size_initial = step_size_initial;
}

// [[Rcpp::export]]
void OdeControl_set_tol_abs(SEXP control_xp, double tol_abs) {
  get_OdeControl(control_xp)->tol_abs = tol_abs;
}

// [[Rcpp::export]]
void OdeControl_set_tol_rel(SEXP control_xp, double tol_rel) {
  get_OdeControl(control_xp)->tol_rel = tol_rel;
}

// [[Rcpp::export]]
void OdeControl_set_a_y(SEXP control_xp, double a_y) {
  get_OdeControl(control_xp)->a_y = a_y;
}

// [[Rcpp::export]]
void OdeControl_set_a_dydt(SEXP control_xp, double a_dydt) {
  get_OdeControl(control_xp)->a_dydt = a_dydt;
}

// [[Rcpp::export]]
void OdeControl_set_step_size_min(SEXP control_xp, double step_size_min) {
  get_OdeControl(control_xp)->step_size_min = step_size_min;
}

// [[Rcpp::export]]
void OdeControl_set_step_size_max(SEXP control_xp, double step_size_max) {
  get_OdeControl(control_xp)->step_size_max = step_size_max;
}

// [[Rcpp::export]]
void OdeControl_set_step_size_initial(SEXP control_xp, double step_size_initial) {
  get_OdeControl(control_xp)->step_size_initial = step_size_initial;
}

//-------------------------------------------------------------------------
// Drivers interface

// [[Rcpp::export]]
SEXP Drivers_new() {
  Rcpp::XPtr<drivers::Drivers> ptr(new drivers::Drivers(), true);
  return ptr;
}

// [[Rcpp::export]]
void Drivers_set_constant(SEXP drivers_xp, std::string driver_name, double k) {
  auto drv = get_Drivers(drivers_xp);
  drv->set_constant(driver_name, k);
}

// [[Rcpp::export]]
void Drivers_set_variable(SEXP drivers_xp, std::string driver_name,
                          Rcpp::NumericVector x, Rcpp::NumericVector y) {
  auto drv = get_Drivers(drivers_xp);
  std::vector<double> x_vec(x.begin(), x.end());
  std::vector<double> y_vec(y.begin(), y.end());
  drv->set_variable(driver_name, x_vec, y_vec);
}

// [[Rcpp::export]]
void Drivers_set_extrapolate(SEXP drivers_xp, std::string driver_name, bool extrapolate) {
  auto drv = get_Drivers(drivers_xp);
  drv->set_extrapolate(driver_name, extrapolate);
}

// [[Rcpp::export]]
double Drivers_evaluate(SEXP drivers_xp, std::string driver_name, double x) {
  auto drv = get_Drivers(drivers_xp);
  return drv->evaluate(driver_name, x);
}

// [[Rcpp::export]]
Rcpp::NumericVector Drivers_evaluate_range(SEXP drivers_xp, std::string driver_name,
                                           Rcpp::NumericVector x) {
  auto drv = get_Drivers(drivers_xp);
  std::vector<double> x_vec(x.begin(), x.end());
  std::vector<double> result = drv->evaluate_range(driver_name, x_vec);
  return Rcpp::wrap(result);
}

// [[Rcpp::export]]
Rcpp::CharacterVector Drivers_get_names(SEXP drivers_xp) {
  auto drv = get_Drivers(drivers_xp);
  std::vector<std::string> names = drv->get_names();
  return Rcpp::wrap(names);
}

// [[Rcpp::export]]
void Drivers_clear(SEXP drivers_xp) {
  auto drv = get_Drivers(drivers_xp);
  drv->clear();
}
