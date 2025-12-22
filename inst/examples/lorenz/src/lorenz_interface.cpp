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
// Assume you already have XPtr<SystemType> and XPtr<ode::OdeControl> in R

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

// Simple accessors
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
  auto y = r->state(); // ode::state_type, typically std::vector<double>
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

// single element as external pointer to a *copy*
// [[Rcpp::export]]
SEXP Solver_get_history_element(SEXP solver_xp, std::size_t i)
{
  auto r = get_solver(solver_xp);
  SystemType elem = r->get_history_element(i); // value copy
  Rcpp::XPtr<SystemType> xp(new SystemType(elem), true);
  return xp;
}

// full history as a list of external pointers
// [[Rcpp::export]]
Rcpp::List Solver_get_history(SEXP solver_xp)
{
  auto r = get_solver(solver_xp);
  std::vector<SystemType> hist = r->get_history();
  R_xlen_t n = static_cast<R_xlen_t>(hist.size());
  Rcpp::List out(n);

  for (R_xlen_t i = 0; i < n; ++i)
  {
    Rcpp::XPtr<SystemType> xp(new SystemType(hist[i]), true);
    out[i] = xp;
  }

  return out;
}

//-------------------------------------------------------------------------
// Rcpp interface for System

// Convenience accessor
inline Rcpp::XPtr<SystemType> get_System(SEXP xp) {
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
  auto lor = get_System(system_xp);
  std::vector<double> p = lor->pars();
  return Rcpp::wrap(p);
}

// State interface

// Set internal state y0, y1, y2 and update rates.
// [[Rcpp::export]]
void System_set_state(SEXP system_xp, Rcpp::NumericVector y) {
  if (y.size() != 3) {
    Rcpp::stop("System_set_state: state vector must have length 3");
  }

  auto lor = get_System(system_xp);

  std::vector<double> tmp(y.begin(), y.end());
  ode::const_iterator it = tmp.begin();
  lor->set_ode_state(it);
}

// Get current state (y0, y1, y2) as numeric(3)
// [[Rcpp::export]]
Rcpp::NumericVector System_state(SEXP system_xp) {
  auto lor = get_System(system_xp);

  std::vector<double> tmp(3);
  ode::iterator it = tmp.begin();
  lor->ode_state(it);

  return Rcpp::wrap(tmp);
}

// Get current rates (dy0dt, dy1dt, dy2dt) as numeric(3)
// [[Rcpp::export]]
Rcpp::NumericVector System_rates(SEXP system_xp) {
  auto lor = get_System(system_xp);

  std::vector<double> tmp(3);
  ode::iterator it = tmp.begin();
  lor->ode_rates(it);

  return Rcpp::wrap(tmp);
}

/*----------------------------------------------------------*/
// Rcpp interface for the ODE control object

// [[Rcpp::export]]
SEXP OdeControl_new(double tol_abs,
                    double tol_rel,
                    double a_y,
                    double a_dydt,
                    double step_size_min,
                    double step_size_max,
                    double step_size_initial)
{
  // allocate new control object
  auto *ctrl = new odelia::ode::OdeControl(
      tol_abs,
      tol_rel,
      a_y,
      a_dydt,
      step_size_min,
      step_size_max,
      step_size_initial);

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
