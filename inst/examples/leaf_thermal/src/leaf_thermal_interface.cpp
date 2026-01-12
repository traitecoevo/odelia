// [[Rcpp::depends(Rcpp, odelia)]]
// [[Rcpp::plugins(cpp20)]]

/* This file defines uses Rcpp to create an interface for the leaf thermal model, including 
- parameter struct: LeafThermalPars
- model system: LeafThermalSystem
- drivers: drivers::Drivers
- ODE solver: ode::Solver<LeafThermalSystem> 

Most are straightforward general wrappers around the underlying C++ methods and not specific to the leaf thermal model.

Throughout the functions below, we 

- use functions names like Class_method to indicate method 'method' for class 'Class'. You should refer to these functions for details on underlying C++ methods.
- use Rcpp::XPtr to manage pointers to C++ objects created in R and passed back to C++.
- convert between Rcpp types (NumericVector, List, DataFrame) and C++ types (std::vector, structs) as needed.
- provide error checking for input sizes and validity where appropriate.
- return results in R-friendly formats.
- all functions are exported to R using the [[Rcpp::export]] attribute.

These functions are intended only as interface and are called via a corresponding R interface,to provide a user-friendly R interface. See `R/leaf_thermal_interface.R` for details.
*/

#include <Rcpp.h>
#include "leaf_thermal_system.hpp"

using namespace Rcpp;
using namespace odelia;

// Define types for convenience
typedef LeafThermalSystem SystemType;
typedef ode::Solver<SystemType> SolverType;

// Helpers to convert pointers from SEXP to XPtr
inline Rcpp::XPtr<SolverType> get_solver(SEXP xp) {
  return Rcpp::XPtr<SolverType>(xp);
}

inline Rcpp::XPtr<LeafThermalSystem> get_LeafThermalSystem(SEXP xp) {
  return Rcpp::XPtr<LeafThermalSystem>(xp);
}

inline Rcpp::XPtr<drivers::Drivers> get_Drivers(SEXP xp) {
  return Rcpp::XPtr<drivers::Drivers>(xp);
}

inline Rcpp::XPtr<ode::OdeControl> get_OdeControl(SEXP xp){
  return Rcpp::XPtr<ode::OdeControl>(xp);
}
//-------------------------------------------------------------------------
// Rcpp interface for ode::Solver<LeafThermalSystem>

// [[Rcpp::export]]
SEXP Solver_new(SEXP system_xp, SEXP control_xp) {
  Rcpp::XPtr<SystemType> sys(system_xp);
  Rcpp::XPtr<ode::OdeControl> ctrl(control_xp);

  Rcpp::XPtr<SolverType> ptr(new SolverType(*sys, *ctrl), true);
  return ptr;
}

// [[Rcpp::export]]
void Solver_reset(SEXP solver_xp) {
  auto r = get_solver(solver_xp);
  r->reset();
}

// [[Rcpp::export]]
double Solver_time(SEXP solver_xp) {
  auto r = get_solver(solver_xp);
  return r->time();
}

// [[Rcpp::export]]
Rcpp::NumericVector Solver_state(SEXP solver_xp) {
  auto r = get_solver(solver_xp);
  auto y = r->state();
  return Rcpp::wrap(y);
}

// [[Rcpp::export]]
Rcpp::NumericVector Solver_times(SEXP solver_xp) {
  auto r = get_solver(solver_xp);
  auto ts = r->times(); // std::vector<double>
  return Rcpp::wrap(ts);
}

// [[Rcpp::export]]
void Solver_set_state(SEXP solver_xp,
                      Rcpp::NumericVector y,
                      double time) {
  auto r = get_solver(solver_xp);
  std::vector<double> yy(y.begin(), y.end());
  r->set_state(yy, time);
}

// [[Rcpp::export]]
void Solver_advance_adaptive(SEXP solver_xp, Rcpp::NumericVector times) {
  auto r = get_solver(solver_xp);
  std::vector<double> ts(times.begin(), times.end());
  r->advance_adaptive(ts);
}

// [[Rcpp::export]]
void Solver_advance_fixed(SEXP solver_xp,
                          Rcpp::NumericVector times) {
  auto r = get_solver(solver_xp);
  std::vector<double> ts(times.begin(), times.end());
  r->advance_fixed(ts);
}

// [[Rcpp::export]]
void Solver_step(SEXP solver_xp) {
  auto r = get_solver(solver_xp);
  r->step();
}

// [[Rcpp::export]]
bool Solver_get_collect(SEXP solver_xp) {
  auto r = get_solver(solver_xp);
  return r->get_collect();
}

// [[Rcpp::export]]
void Solver_set_collect(SEXP solver_xp, bool x) {
  auto r = get_solver(solver_xp);
  r->set_collect(x);
}

// [[Rcpp::export]]
std::size_t Solver_get_history_size(SEXP solver_xp) {
  auto r = get_solver(solver_xp);
  return r->get_history_size();
}

// Return history for a single step as a dataframe
// [[Rcpp::export]]
Rcpp::DataFrame Solver_get_history_step(SEXP solver_xp, std::size_t i) {
  auto r = get_solver(solver_xp);
  int nrows = r->get_history_size();
  if (i >= nrows)
  {
    Rcpp::stop("Index out of bounds. History size is " +
               std::to_string(nrows) + ", requested index " +
               std::to_string(i));
  }

  SystemType elem = r->get_history_step(i);
  std::vector<double> out = elem.record_step();
  std::vector<std::string> names = r->get_system().record_colnames();

  Rcpp::List df_list(names.size());
  for (size_t j = 0; j < names.size(); ++j)
  {
    df_list[j] = Rcpp::wrap(out[j]);
  }
  df_list.attr("names") = names;

  return Rcpp::DataFrame(df_list);
}

// Return full history as a dataframe
// [[Rcpp::export]]
Rcpp::List Solver_get_history(SEXP solver_xp) {
  auto r = get_solver(solver_xp);
  int nrows = r->get_history_size();

  std::vector<std::string> names = r->get_system().record_colnames();
  int ncols = names.size();

  // Prepare column storage and reserve rows if you can estimate nrows
  std::vector<std::vector<double>> cols(ncols);
  for (size_t j = 0; j < ncols; ++j)
  {
    cols[j].reserve(nrows);
  }

  // Produce rows and push into columns
  for (size_t i = 0; i < nrows; ++i)
  {
    SystemType elem = r->get_history_step(i);
    std::vector<double> row = elem.record_step();

    for (size_t j = 0; j < ncols; ++j)
    {
      cols[j].push_back(row[j]);
    }
  }

  // Convert each std::vector<double> -> Rcpp::NumericVector
  // Construct an R list with named columns, then DataFrame::create from it
  Rcpp::List out(ncols);

  for (size_t j = 0; j < ncols; ++j)
  {
    // Efficient conversion: construct NumericVector from iterators (copies once)
    NumericVector nv(cols[j].begin(), cols[j].end());
    out[j] = nv;
  }

  // Set column names
  out.attr("names") = names;

  // Return a DataFrame
  // NB (DataFrame::create duplicates columns again if you pass names differently; constructing directly from a named List is more efficient)

  return DataFrame(out);
}

//-------------------------------------------------------------------------
// Rcpp interface for LeafThermalSystem

// helper: C++ struct -> R list
inline List leaf_pars_to_list(const LeafThermalPars &p) {
  return List::create(
      _["k_H"] = p.k_H,
      _["g_tr_max"] = p.g_tr_max,
      _["m_tr"] = p.m_tr,
      _["T_tr_mid"] = p.T_tr_mid);
}

// helper: R list -> C++ struct (missing entries -> keep default)
inline LeafThermalPars leaf_pars_from_list(const List &L) {
  LeafThermalPars p; // starts with C++ defaults

  if (L.containsElementNamed("k_H"))
    p.k_H = as<double>(L["k_H"]);
  if (L.containsElementNamed("g_tr_max"))
    p.g_tr_max = as<double>(L["g_tr_max"]);
  if (L.containsElementNamed("m_tr"))
    p.m_tr = as<double>(L["m_tr"]);
  if (L.containsElementNamed("T_tr_mid"))
    p.T_tr_mid = as<double>(L["T_tr_mid"]);

  return p;
}

// [[Rcpp::export]]
Rcpp::List LeafThermalSystemPars() {
  LeafThermalPars p; // uses C++ defaults
  return leaf_pars_to_list(p);
}

// [[Rcpp::export]]
SEXP LeafThermalSystem_new(Rcpp::List pars_list, SEXP drivers_xp) {

  LeafThermalPars pars = leaf_pars_from_list(pars_list);
  auto drivers = get_Drivers(drivers_xp);
  // R will delete this when the external pointer is GC'd
  Rcpp::XPtr<LeafThermalSystem> ptr(new LeafThermalSystem(pars, *drivers), true);
  return ptr;
}

// [[Rcpp::export]]
Rcpp::NumericVector LeafThermalSystem_pars(SEXP LeafThermalSystem_xp) {
  auto lor = get_LeafThermalSystem(LeafThermalSystem_xp);
  std::vector<double> p = lor->get_pars();
  return Rcpp::wrap(p);
}

// [[Rcpp::export]]
void LeafThermalSystem_initialize_drivers(SEXP LeafThermalSystem_xp, SEXP drivers_xp) {

  auto drv = get_Drivers(drivers_xp);
  auto lor = get_LeafThermalSystem(LeafThermalSystem_xp);
  
  lor->initialize_drivers(*drv);
}

// [[Rcpp::export]]
void LeafThermalSystem_set_state(SEXP LeafThermalSystem_xp, Rcpp::NumericVector y, double time) {

  auto lor = get_LeafThermalSystem(LeafThermalSystem_xp);

  if (y.size() != lor->ode_size())
  {
    Rcpp::stop("LeafThermalSystem_set_state: state vector must have length 1");
  }

  std::vector<double> tmp(y.begin(), y.end());
  ode::const_iterator it = tmp.begin();
  lor->set_ode_state(it, time);
}

// [[Rcpp::export]]
Rcpp::NumericVector LeafThermalSystem_state(SEXP LeafThermalSystem_xp) {
  auto lor = get_LeafThermalSystem(LeafThermalSystem_xp);

  std::vector<double> tmp(lor->ode_size());
  ode::iterator it = tmp.begin();
  lor->ode_state(it);

  return Rcpp::wrap(tmp);
}

// [[Rcpp::export]]
Rcpp::NumericVector LeafThermalSystem_rates(SEXP LeafThermalSystem_xp) {
  auto lor = get_LeafThermalSystem(LeafThermalSystem_xp);

  // vectors to hold rates
  std::vector<double> tmp(lor->ode_size());
  // fill rates using iterator
  ode::iterator it = tmp.begin();
  lor->ode_rates(it);

  // package and return
  return Rcpp::wrap(tmp);
}

// [[Rcpp::export]]
Rcpp::NumericVector LeafThermalSystem_get_current_drivers(SEXP LeafThermalSystem_xp) {
  auto lor = get_LeafThermalSystem(LeafThermalSystem_xp);
  return Rcpp::wrap(lor->get_current_drivers());
}

//-------------------------------------------------------------------------
// Rcpp interface for odelia::ode::OdeControl

// [[Rcpp::export]]
SEXP OdeControl_new()
{
  Rcpp::XPtr<ode::OdeControl> ctrl(new ode::OdeControl(), true);
  return ctrl;
}

// [[Rcpp::export]]
void OdeControl_set_controls(SEXP control_xp, double atol, double rtol, 
          double a, double b, double eps, double h, double hmax) {
  auto ctrl = get_OdeControl(control_xp);
  ctrl->set_controls(atol, rtol, a, b, eps, h, hmax);
}

// [[Rcpp::export]]
Rcpp::NumericVector OdeControl_get_controls(SEXP control_xp) {
  auto ctrl = get_OdeControl(control_xp);
  return Rcpp::wrap(ctrl->get_controls());
}

// [[Rcpp::export]]
void OdeControl_set_tol_abs(SEXP control_xp, double atol) {
  auto ctrl = get_OdeControl(control_xp);
  ctrl->set_tol_abs(atol);
}

// [[Rcpp::export]]
void OdeControl_set_tol_rel(SEXP control_xp, double rtol) {
  auto ctrl = get_OdeControl(control_xp);
  ctrl->set_tol_rel(rtol);
}

// [[Rcpp::export]]
void OdeControl_set_a_y(SEXP control_xp, double a_y) {
  auto ctrl = get_OdeControl(control_xp);
  ctrl->set_a_y(a_y);
}

// [[Rcpp::export]]
void OdeControl_set_a_dydt(SEXP control_xp, double a_dydt) {
  auto ctrl = get_OdeControl(control_xp);
  ctrl->set_a_dydt(a_dydt);
}

// [[Rcpp::export]]
void OdeControl_set_step_size_min(SEXP control_xp, double step_size_min) {
  auto ctrl = get_OdeControl(control_xp);
  ctrl->set_step_size_min(step_size_min);
}

// [[Rcpp::export]]
void OdeControl_set_step_size_max(SEXP control_xp, double step_size_max) {
  auto ctrl = get_OdeControl(control_xp);
  ctrl->set_step_size_max(step_size_max);
}

// [[Rcpp::export]]
void OdeControl_set_step_size_initial(SEXP control_xp, double step_size_initial) {
  auto ctrl = get_OdeControl(control_xp);
  ctrl->set_step_size_initial(step_size_initial);
}

//-------------------------------------------------------------------------
// Rcpp interface for Drivers class

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
