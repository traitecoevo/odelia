/* Generic ODE interface for odelia package
 * 
 * This file provides generic functions that can be used with any ODE system.
 * It includes:
 * - Drivers: External forcing data (constant or time-varying)
 * - OdeControl: Solver settings and tolerances
 * 
 * These functions are system-agnostic and can be used by any model
 * (Lorenz, Leaf Thermal, ATLAS, etc.)
 */

#include <Rcpp.h>
#include <odelia/drivers.hpp>
#include <odelia/ode_control.hpp>

using namespace Rcpp;
using namespace odelia;

//-------------------------------------------------------------------------
// OdeControl interface

// [[Rcpp::export]]
SEXP OdeControl_new() {
  return Rcpp::XPtr<ode::OdeControl>(new ode::OdeControl(), true);
}

//-------------------------------------------------------------------------
// Drivers interface

// Helper to get Drivers pointer
inline Rcpp::XPtr<drivers::Drivers> get_Drivers(SEXP xp) {
  return Rcpp::XPtr<drivers::Drivers>(xp);
}

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