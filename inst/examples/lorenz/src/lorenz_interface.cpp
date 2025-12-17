// [[Rcpp::depends(Rcpp, odelia)]]
// [[Rcpp::plugins(cpp20)]]

#include <Rcpp.h>
#include "lorenz_solver.hpp"

using namespace Rcpp;
using namespace odelia;

// [[Rcpp::export]]
std::vector<double> get_number()
{
  ode::LorenzSystem sys{10.0, 28.0, 8.0 / 3.0};

  ode::OdeControl ctrl(1e-8, 1e-8, 1.0, 0.0,
                  1e-8, 10.0, 1e-6);

  std::cout <<sys.ode_size() << std::endl;

  ode::Solver<ode::LorenzSystem> solver(sys, ctrl);
  std::cout << solver.get_time() << std::endl;

  solver.advance_adaptive(sys, 10000);
  std::vector<double> times = solver.get_times();

  return times;
}

//-------------------------------------------------------------------------
// Rcpp interface for the ODE Runner 

// Helper to get the XPtr<Runner>
inline Rcpp::XPtr<ode::Runner> get_runner(SEXP xp) {
  return Rcpp::XPtr<ode::Runner>(xp);
}

// Constructor: build Runner from existing LorenzSystem + OdeControl 
// Assume you already have XPtr<LorenzSystem> and XPtr<ode::OdeControl> in R

// [[Rcpp::export]]
SEXP Runner_new(SEXP LorenzSystem_xp, SEXP control_xp)
{
  Rcpp::XPtr<ode::LorenzSystem> sys(LorenzSystem_xp);
  Rcpp::XPtr<ode::OdeControl> ctrl(control_xp);

  Rcpp::XPtr<ode::Runner> ptr(new ode::Runner(*sys, *ctrl), true);
  return ptr;
}

// [[Rcpp::export]]
void Runner_reset(SEXP runner_xp)
{
  auto r = get_runner(runner_xp);
  r->reset();
}

// Simple accessors
// [[Rcpp::export]]
double Runner_time(SEXP runner_xp)
{
  auto r = get_runner(runner_xp);
  return r->time();
}

// [[Rcpp::export]]
Rcpp::NumericVector Runner_state(SEXP runner_xp)
{
  auto r = get_runner(runner_xp);
  auto y = r->state(); // ode::state_type, typically std::vector<double>
  return Rcpp::wrap(y);
}

// [[Rcpp::export]]
Rcpp::NumericVector Runner_times(SEXP runner_xp)
{
  auto r = get_runner(runner_xp);
  auto ts = r->times(); // std::vector<double>
  return Rcpp::wrap(ts);
}

// 3. Mutators / evolution
// [[Rcpp::export]]
void Runner_set_state(SEXP runner_xp,
                      Rcpp::NumericVector y,
                      double time)
{
  auto r = get_runner(runner_xp);
  std::vector<double> yy(y.begin(), y.end());
  r->set_state(yy, time);
}


// [[Rcpp::export]]
void Runner_advance_adaptive(SEXP runner_xp, Rcpp::NumericVector times)
{
  auto r = get_runner(runner_xp);
  std::vector<double> ts(times.begin(), times.end());
  r->advance_adaptive(ts);
}

// [[Rcpp::export]]
void Runner_advance_fixed(SEXP runner_xp,
                          Rcpp::NumericVector times)
{
  auto r = get_runner(runner_xp);
  std::vector<double> ts(times.begin(), times.end());
  r->advance_fixed(ts);
}

// [[Rcpp::export]]
void Runner_step(SEXP runner_xp)
{
  auto r = get_runner(runner_xp);
  r->step();
}

//-------------------------
// collect flag
//-------------------------

// [[Rcpp::export]]
bool Runner_get_collect(SEXP runner_xp)
{
  auto r = get_runner(runner_xp);
  return r->get_collect();
}

// [[Rcpp::export]]
void Runner_set_collect(SEXP runner_xp, bool x)
{
  auto r = get_runner(runner_xp);
  r->set_collect(x);
}

//-------------------------
// history access
//-------------------------

// size of history
// [[Rcpp::export]]
std::size_t Runner_get_history_size(SEXP runner_xp)
{
  auto r = get_runner(runner_xp);
  return r->get_history_size();
}

// single element as external pointer to a *copy*
// [[Rcpp::export]]
SEXP Runner_get_history_element(SEXP runner_xp, std::size_t i)
{
  auto r = get_runner(runner_xp);
  ode::LorenzSystem elem = r->get_history_element(i); // value copy
  Rcpp::XPtr<ode::LorenzSystem> xp(new ode::LorenzSystem(elem), true);
  return xp;
}

// full history as a list of external pointers
// [[Rcpp::export]]
Rcpp::List Runner_get_history(SEXP runner_xp)
{
  auto r = get_runner(runner_xp);
  std::vector<ode::LorenzSystem> hist = r->get_history();

  R_xlen_t n = static_cast<R_xlen_t>(hist.size());
  Rcpp::List out(n);

  for (R_xlen_t i = 0; i < n; ++i)
  {
    Rcpp::XPtr<ode::LorenzSystem> xp(new ode::LorenzSystem(hist[i]), true);
    out[i] = xp;
  }

  return out;
}

//-------------------------------------------------------------------------
// Rcpp interface for LorenzSystem

// Convenience accessor
inline Rcpp::XPtr<ode::LorenzSystem> get_LorenzSystem(SEXP xp) {
  return Rcpp::XPtr<ode::LorenzSystem>(xp);
}

// Constructors  & basic access
// [[Rcpp::export]]
SEXP LorenzSystem_new(double sigma, double R, double b) {
  // R will delete this when the external pointer is GC'd
  Rcpp::XPtr<ode::LorenzSystem> ptr(new ode::LorenzSystem(sigma, R, b), true);
  return ptr;
}

// [[Rcpp::export]]
Rcpp::NumericVector LorenzSystem_pars(SEXP LorenzSystem_xp) {
  auto lor = get_LorenzSystem(LorenzSystem_xp);
  std::vector<double> p = lor->pars();
  return Rcpp::wrap(p);
}

// State interface

// Set internal state y0, y1, y2 and update rates.
// [[Rcpp::export]]
void LorenzSystem_set_state(SEXP LorenzSystem_xp, Rcpp::NumericVector y) {
  if (y.size() != 3) {
    Rcpp::stop("LorenzSystem_set_state: state vector must have length 3");
  }

  auto lor = get_LorenzSystem(LorenzSystem_xp);

  std::vector<double> tmp(y.begin(), y.end());
  ode::const_iterator it = tmp.begin();
  lor->set_ode_state(it);
}

// Get current state (y0, y1, y2) as numeric(3)
// [[Rcpp::export]]
Rcpp::NumericVector LorenzSystem_state(SEXP LorenzSystem_xp) {
  auto lor = get_LorenzSystem(LorenzSystem_xp);

  std::vector<double> tmp(3);
  ode::iterator it = tmp.begin();
  lor->ode_state(it);

  return Rcpp::wrap(tmp);
}

// Get current rates (dy0dt, dy1dt, dy2dt) as numeric(3)
// [[Rcpp::export]]
Rcpp::NumericVector LorenzSystem_rates(SEXP LorenzSystem_xp) {
  auto lor = get_LorenzSystem(LorenzSystem_xp);

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
// Define single function that gives rates of change for the Lorenz system 
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
