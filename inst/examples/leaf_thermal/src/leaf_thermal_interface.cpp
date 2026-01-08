// [[Rcpp::depends(Rcpp, odelia)]]
// [[Rcpp::plugins(cpp20)]]

#include <Rcpp.h>
#include <odelia/ode_solver.hpp>
#include "leaf_thermal_system.hpp"

using namespace Rcpp;
using namespace odelia;
using odelia::ode::LeafThermalPars;

// Rcpp interface for the ODE Solver
typedef odelia::ode::LeafThermalSystem SystemType;
typedef ode::Solver<SystemType> SolverType;

//-------------------------------------------------------------------------
// Rcpp interface for the ODE Runner


// Helper to get the XPtr<Solver>
inline Rcpp::XPtr<SolverType> get_solver(SEXP xp)
{
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
  std::vector<std::string> names = r->get_system().record_colnames();

  Rcpp::List df_list(names.size());
  for (size_t j = 0; j < names.size(); ++j)
  {
    df_list[j] = Rcpp::wrap(out[j]);
  }
  df_list.attr("names") = names;

  return Rcpp::DataFrame(df_list);
}

// full history as a dataframe
// [[Rcpp::export]]
Rcpp::List Solver_get_history(SEXP solver_xp)
{
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
inline List leaf_pars_to_list(const LeafThermalPars &p)
{
  return List::create(
      _["Tmean"] = p.Tmean,
      _["Tamp"] = p.Tamp,
      _["tpeak"] = p.tpeak,
      _["k_H"] = p.k_H,
      _["g_tr_max"] = p.g_tr_max,
      _["m_tr"] = p.m_tr,
      _["T_tr_mid"] = p.T_tr_mid);
}

// helper: R list -> C++ struct (missing entries -> keep default)
inline LeafThermalPars leaf_pars_from_list(const List &L)
{
  LeafThermalPars p; // starts with C++ defaults

  if (L.containsElementNamed("Tmean"))
    p.Tmean = as<double>(L["Tmean"]);
  if (L.containsElementNamed("Tamp"))
    p.Tamp = as<double>(L["Tamp"]);
  if (L.containsElementNamed("tpeak"))
    p.tpeak = as<double>(L["tpeak"]);

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
Rcpp::List leaf_thermal_pars()
{
  LeafThermalPars p; // uses C++ defaults
  return leaf_pars_to_list(p);
}

// [[Rcpp::export]]
Rcpp::List leaf_thermal_pars_complete(Rcpp::List pars)
{
  // Fill missing entries with defaults, return completed list
  LeafThermalPars p = leaf_pars_from_list(pars);
  return leaf_pars_to_list(p);
}

// Convenience accessor
inline Rcpp::XPtr<ode::LeafThermalSystem> get_LeafThermalSystem(SEXP xp) {
  return Rcpp::XPtr<ode::LeafThermalSystem>(xp);
}

// Constructors  & basic access
// [[Rcpp::export]]
SEXP LeafThermalSystem_new(Rcpp::List pars_list) {

  LeafThermalPars pars = leaf_pars_from_list(pars_list);
  // R will delete this when the external pointer is GC'd
  Rcpp::XPtr<ode::LeafThermalSystem> ptr(new ode::LeafThermalSystem(pars), true);
  return ptr;
}


// [[Rcpp::export]]
Rcpp::NumericVector LeafThermalSystem_pars(SEXP LeafThermalSystem_xp) {
  auto lor = get_LeafThermalSystem(LeafThermalSystem_xp);
  std::vector<double> p = lor->get_pars();
  return Rcpp::wrap(p);
}



// State interface

// Set internal state and update rates.
// [[Rcpp::export]]
void LeafThermalSystem_set_state(SEXP LeafThermalSystem_xp, Rcpp::NumericVector y) {
  if (y.size() != 1) {
    Rcpp::stop("LeafThermalSystem_set_state: state vector must have length 1");
  }

  auto lor = get_LeafThermalSystem(LeafThermalSystem_xp);

  std::vector<double> tmp(y.begin(), y.end());
  ode::const_iterator it = tmp.begin();
  lor->set_ode_state(it);
}

// Get current state (y0, y1, y2) as numeric(3)
// [[Rcpp::export]]
Rcpp::NumericVector LeafThermalSystem_state(SEXP LeafThermalSystem_xp) {
  auto lor = get_LeafThermalSystem(LeafThermalSystem_xp);

  std::vector<double> tmp(1);
  ode::iterator it = tmp.begin();
  lor->ode_state(it);

  return Rcpp::wrap(tmp);
}

// [[Rcpp::export]]
Rcpp::NumericVector LeafThermalSystem_rates(SEXP LeafThermalSystem_xp) {
  auto lor = get_LeafThermalSystem(LeafThermalSystem_xp);

  std::vector<double> tmp(1);
  ode::iterator it = tmp.begin();
  lor->ode_rates(it);

  return Rcpp::wrap(tmp);
}

/*----------------------------------------------------------*/
// Rcpp interface for the ODE control object

// [[Rcpp::export]]
SEXP OdeControl_default()
{
  // allocate new control object
  auto *ctrl = new odelia::ode::OdeControl(1e-8, 1e-8, 1.0, 1.0, 1e-8, 1.0, 1e-3);

  // wrap in an external pointer, auto-delete on GC
  Rcpp::XPtr<odelia::ode::OdeControl> xp(ctrl, true);
  return xp;
}
