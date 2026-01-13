/* This file defines uses Rcpp to create an interface for the lorenz model, including

- model system: LorenzSystem
- ODE solver: ode::Solver<LorenzSystem>

Most are straightforward general wrappers around the underlying C++ methods and not specific to the lorenz model.

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
#include <odelia/ode_solver.hpp>
#include <examples/lorenz_system.hpp>

using namespace Rcpp;
using namespace odelia;

// Define types for convenience
typedef LorenzSystem SystemType;
typedef ode::Solver<SystemType> SolverType;

// Helpers to convert pointers from SEXP to XPtr
inline Rcpp::XPtr<SolverType> get_solver(SEXP xp) {
  return Rcpp::XPtr<SolverType>(xp);
}

inline Rcpp::XPtr<SystemType> get_system(SEXP xp) {
  return Rcpp::XPtr<SystemType>(xp);
}

//-------------------------------------------------------------------------
// Rcpp interface for ode::Solver<LorenzSystem>

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

// Helper to get column names
CharacterVector get_column_names() {
  CharacterVector names(7);
  names = {"time", "x", "y", "z", "dxdt", "dydt", "dzdt"};
  return names;
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
  
  CharacterVector names = get_column_names();

  Rcpp::List df_list(names.size());
  for (size_t j = 0; j < names.size(); ++j) {
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

  CharacterVector names = get_column_names();
  int ncols = names.size();

  //Prepare column storage and reserve rows if you can estimate nrows
  std::vector<std::vector<double>> cols(ncols);
  for (size_t j = 0; j < ncols; ++j) {
    cols[j].reserve(nrows);
  }

  // Produce rows and push into columns
  for (size_t i = 0; i < nrows; ++i) {
    SystemType elem = r->get_history_step(i);
    std::vector<double> row = elem.record_step();

    for (size_t j = 0; j < ncols; ++j) {
      cols[j].push_back(row[j]);
    }
  }

  // Convert each std::vector<double> -> Rcpp::NumericVector
  // Construct an R list with named columns, then DataFrame::create from it
  Rcpp::List out(ncols);
  
  for (size_t j = 0; j < ncols; ++j) {
    // Efficient conversion: construct NumericVector from iterators (copies once)
    NumericVector nv(cols[j].begin(), cols[j].end());
    out[j] = nv;
  }
  
  // Set column names
  out.attr("names") = get_column_names();

  // Return a DataFrame 
  // NB (DataFrame::create duplicates columns again if you pass names differently; constructing directly from a named List is more efficient)
  
  return DataFrame(out);
}

//-------------------------------------------------------------------------
// Rcpp interface for System

// Constructors  & basic access
// [[Rcpp::export]]
SEXP System_new(double sigma, double R, double b) {
  // R will delete this when the external pointer is GC'd
  Rcpp::XPtr<SystemType> ptr(new SystemType(sigma, R, b), true);
  return ptr;
}

// [[Rcpp::export]]
Rcpp::NumericVector System_pars(SEXP system_xp) {
  auto lor = get_system(system_xp);
  std::vector<double> p = lor->pars();
  return Rcpp::wrap(p);
}

// [[Rcpp::export]]
void System_set_state(SEXP system_xp, Rcpp::NumericVector y, double time) {
  
  auto lor = get_system(system_xp);

  size_t ode_size = lor->ode_size();
  if (y.size() != ode_size)
  {
    Rcpp::stop("System_set_state: state vector must have length " + std::to_string(ode_size));
  }

  std::vector<double> tmp(y.begin(), y.end());
  ode::const_iterator it = tmp.begin();
  lor->set_ode_state(it, time);
}

// [[Rcpp::export]]
Rcpp::NumericVector System_state(SEXP system_xp) {
  auto lor = get_system(system_xp);

  size_t ode_size = lor->ode_size();
  std::vector<double> tmp(ode_size);
  ode::iterator it = tmp.begin();
  lor->ode_state(it);

  return Rcpp::wrap(tmp);
}

// [[Rcpp::export]]
Rcpp::NumericVector System_rates(SEXP system_xp) {
  auto lor = get_system(system_xp);

  size_t ode_size = lor->ode_size();
  std::vector<double> tmp(ode_size);

  ode::iterator it = tmp.begin();
  lor->ode_rates(it);

  return Rcpp::wrap(tmp);
}

//----------------------------------------------------------
// Rcpp interface for the ODE control object

// [[Rcpp::export]]
SEXP OdeControl_new() {
  // allocate new control object
  auto *ctrl = new odelia::ode::OdeControl();

  // wrap in an external pointer, auto-delete on GC
  Rcpp::XPtr<odelia::ode::OdeControl> xp(ctrl, true);
  return xp;
}

//-------------------------------------------------------------------------
// For comparison: Define single function that gives rates of change for the  system 
// that could be passed into deSolve, executed through C++ for speed.

// [[Rcpp::export]]
List lorenz_rhs(double t, NumericVector state, NumericVector pars) {
  double x = state["x"];
  double y = state["y"];
  double z = state["z"];

  double sigma = pars["sigma"];
  double rho = pars["rho"];
  double beta = pars["beta"];

  double dx = sigma * (y - x);
  double dy = x * (rho - z) - y;
  double dz = x * y - beta * z;

  return List::create(NumericVector::create(dx, dy, dz));
}
