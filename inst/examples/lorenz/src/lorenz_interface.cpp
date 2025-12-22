// [[Rcpp::depends(Rcpp, odelia)]]
// [[Rcpp::plugins(cpp20)]]

#include <Rcpp.h>
#include <odelia/ode_solver.hpp>
#include "lorenz_system.hpp"

using namespace Rcpp;
using namespace odelia;

// Rcpp interface for the ODE Solver
typedef odelia::ode::LorenzSystem SystemType;
typedef ode::Solver<SystemType> SolverType;

// Helper to get the XPtr<Solver>
inline Rcpp::XPtr<SolverType> get_solver(SEXP xp) {
  return Rcpp::XPtr<SolverType>(xp);
}

// Constructor: build Solver from existing System + OdeControl
// Assume you already have these defined and pass in as external pointers

// [[Rcpp::export]]
SEXP Solver_new(SEXP system_xp, SEXP control_xp)
{
  Rcpp::XPtr<SystemType> sys(system_xp);
  Rcpp::XPtr<ode::OdeControl> ctrl(control_xp);

  Rcpp::XPtr<SolverType> ptr(new SolverType(*sys, *ctrl), true);
  return ptr;
}

// [[Rcpp::export]]
void Solver_reset(SEXP solver_xp)
{
  auto r = get_solver(solver_xp);
  r->reset();
}

// [[Rcpp::export]]
double Solver_time(SEXP solver_xp)
{
  auto r = get_solver(solver_xp);
  return r->time();
}

// [[Rcpp::export]]
Rcpp::NumericVector Solver_state(SEXP solver_xp)
{
  auto r = get_solver(solver_xp);
  auto y = r->state();
  return Rcpp::wrap(y);
}

// [[Rcpp::export]]
Rcpp::NumericVector Solver_times(SEXP solver_xp)
{
  auto r = get_solver(solver_xp);
  auto ts = r->times(); // std::vector<double>
  return Rcpp::wrap(ts);
}

// 3. Mutators / evolution
// [[Rcpp::export]]
void Solver_set_state(SEXP solver_xp,
                      Rcpp::NumericVector y,
                      double time)
{
  auto r = get_solver(solver_xp);
  std::vector<double> yy(y.begin(), y.end());
  r->set_state(yy, time);
}

// [[Rcpp::export]]
void Solver_advance_adaptive(SEXP solver_xp, Rcpp::NumericVector times)
{
  auto r = get_solver(solver_xp);
  std::vector<double> ts(times.begin(), times.end());
  r->advance_adaptive(ts);
}

// [[Rcpp::export]]
void Solver_advance_fixed(SEXP solver_xp,
                          Rcpp::NumericVector times)
{
  auto r = get_solver(solver_xp);
  std::vector<double> ts(times.begin(), times.end());
  r->advance_fixed(ts);
}

// [[Rcpp::export]]
void Solver_step(SEXP solver_xp)
{
  auto r = get_solver(solver_xp);
  r->step();
}

//-------------------------
// collect flag
//-------------------------

// [[Rcpp::export]]
bool Solver_get_collect(SEXP solver_xp)
{
  auto r = get_solver(solver_xp);
  return r->get_collect();
}

// [[Rcpp::export]]
void Solver_set_collect(SEXP solver_xp, bool x)
{
  auto r = get_solver(solver_xp);
  r->set_collect(x);
}

//-------------------------
// history access
//-------------------------

// size of history
// [[Rcpp::export]]
std::size_t Solver_get_history_size(SEXP solver_xp)
{
  auto r = get_solver(solver_xp);
  return r->get_history_size();
}

CharacterVector get_column_names() {
  CharacterVector names(7);
  names = {"time", "x", "y", "z", "dxdt", "dydt", "dzdt"};
  return names;
}

// single step
// [[Rcpp::export]]
Rcpp::DataFrame Solver_get_history_step(SEXP solver_xp, std::size_t i)
{
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

// full history as a list of external pointers
// [[Rcpp::export]]
Rcpp::List Solver_get_history(SEXP solver_xp)
{
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

  // 4) Return a DataFrame 
  // NB (DataFrame::create duplicates columns again if you pass names differently; constructing directly from a named List is more efficient)
  
  return DataFrame(out);
}

//-------------------------------------------------------------------------
// Rcpp interface for System

// Convenience accessor
inline Rcpp::XPtr<SystemType> get_system(SEXP xp) {
  return Rcpp::XPtr<SystemType>(xp);
}

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

// State interface

// Set internal state and update rates.
// [[Rcpp::export]]
void System_set_state(SEXP system_xp, Rcpp::NumericVector y) {
  
  auto lor = get_system(system_xp);

  size_t ode_size = lor->ode_size();
  if (y.size() != ode_size)
  {
    Rcpp::stop("System_set_state: state vector must have length " + std::to_string(ode_size));
  }

  std::vector<double> tmp(y.begin(), y.end());
  ode::const_iterator it = tmp.begin();
  lor->set_ode_state(it);
}

// Get current state
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

/*----------------------------------------------------------*/
// Rcpp interface for the ODE control object

// [[Rcpp::export]]
SEXP OdeControl_new()
{
  // allocate new control object
  auto *ctrl = new odelia::ode::OdeControl();

  // wrap in an external pointer, auto-delete on GC
  Rcpp::XPtr<odelia::ode::OdeControl> xp(ctrl, true);
  return xp;
}

//-------------------------------------------------------------------------
// Define single function that gives rates of change for the  system 
// that could be passed into deSolve, executed through C++ for speed.

// [[Rcpp::export]]
List lorenz_rhs(double t, NumericVector state, NumericVector pars)
{
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
