// [[Rcpp::depends(Rcpp, odestepper)]]
// [[Rcpp::plugins(cpp20)]]

#include <Rcpp.h>
#include "../include/lorenz_solver.hpp"

using namespace Rcpp;
using namespace odestepper;

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

// [[Rcpp::export]]
std::vector<double> get_number()
{
  ode::Lorenz sys{10.0, 28.0, 8.0 / 3.0};

  ode::OdeControl ctrl(1e-8, 1e-8, 1.0, 0.0,
                  1e-8, 10.0, 1e-6);

  std::cout <<sys.ode_size() << std::endl;

  ode::Solver<ode::Lorenz> solver(sys, ctrl);

  std::cout << solver.get_time() << std::endl;

  solver.advance_adaptive(sys, 10000);
  std::vector<double> times = solver.get_times();

  std::cout << solver.get_time() << " " << times.size() << std::endl;

  return times;
}

// Helper to get the XPtr<Runner> ---------------------------------------------

inline Rcpp::XPtr<ode::Runner> get_runner(SEXP xp)
{
  return Rcpp::XPtr<ode::Runner>(xp);
}

// 1. Constructor: build Runner from existing Lorenz + OdeControl  ------------
// Here we assume you already have XPtr<Lorenz> and XPtr<ode::OdeControl> in R

// [[Rcpp::export]]
SEXP Runner_new(SEXP lorenz_xp, SEXP control_xp)
{
  Rcpp::XPtr<ode::Lorenz> obj(lorenz_xp);
  Rcpp::XPtr<ode::OdeControl> ctrl(control_xp);

  Rcpp::XPtr<ode::Runner> ptr(new ode::Runner(*obj, *ctrl), true);
  return ptr;
}

// 2. Simple accessors --------------------------------------------------------

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

// 3. Mutators / evolution ----------------------------------------------------

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
void Runner_set_state_from_system(SEXP runner_xp)
{
  auto r = get_runner(runner_xp);
  r->set_state_from_system();
}

// [[Rcpp::export]]
void Runner_advance_adaptive(SEXP runner_xp, double time)
{
  auto r = get_runner(runner_xp);
  r->advance_adaptive(time);
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

// [[Rcpp::export]]
void Runner_step_to(SEXP runner_xp, double time)
{
  auto r = get_runner(runner_xp);
  r->step_to(time);
}

// Convenience accessor
inline Rcpp::XPtr<ode::Lorenz> get_lorenz(SEXP xp) {
  return Rcpp::XPtr<ode::Lorenz>(xp);
}

//-----------------------------
// Constructors / basic access
//-----------------------------

// [[Rcpp::export]]
SEXP Lorenz_new(double sigma, double R, double b) {
  // R will delete this when the external pointer is GC'd
  Rcpp::XPtr<ode::Lorenz> ptr(new ode::Lorenz(sigma, R, b), true);
  return ptr;
}

// [[Rcpp::export]]
Rcpp::NumericVector Lorenz_pars(SEXP lorenz_xp) {
  auto lor = get_lorenz(lorenz_xp);
  std::vector<double> p = lor->pars();
  return Rcpp::wrap(p);
}

//-----------------------------
// State interface
//-----------------------------

// Set internal state y0, y1, y2 and update rates.
// y must be length 3.
// [[Rcpp::export]]
void Lorenz_set_state(SEXP lorenz_xp, Rcpp::NumericVector y) {
  if (y.size() != 3) {
    Rcpp::stop("Lorenz_set_state: state vector must have length 3");
  }

  auto lor = get_lorenz(lorenz_xp);

  std::vector<double> tmp(y.begin(), y.end());
  ode::const_iterator it = tmp.begin();
  lor->set_ode_state(it);
}

// Get current state (y0, y1, y2) as numeric(3)
// [[Rcpp::export]]
Rcpp::NumericVector Lorenz_state(SEXP lorenz_xp) {
  auto lor = get_lorenz(lorenz_xp);

  std::vector<double> tmp(3);
  ode::iterator it = tmp.begin();
  lor->ode_state(it);

  return Rcpp::wrap(tmp);
}

// Get current rates (dy0dt, dy1dt, dy2dt) as numeric(3)
// [[Rcpp::export]]
Rcpp::NumericVector Lorenz_rates(SEXP lorenz_xp) {
  auto lor = get_lorenz(lorenz_xp);

  std::vector<double> tmp(3);
  ode::iterator it = tmp.begin();
  lor->ode_rates(it);

  return Rcpp::wrap(tmp);
}

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
  auto *ctrl = new odestepper::ode::OdeControl(
      tol_abs,
      tol_rel,
      a_y,
      a_dydt,
      step_size_min,
      step_size_max,
      step_size_initial);

  // wrap in an external pointer, auto-delete on GC
  Rcpp::XPtr<odestepper::ode::OdeControl> xp(ctrl, true);
  return xp;
}
