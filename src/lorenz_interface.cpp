/* Lorenz system interface for odelia package
 * 
 * This file provides Lorenz-specific functions:
 * - LorenzSystem creation and manipulation
 * - Solver creation (calls generic templates with LorenzSystem type)
 * 
 * Most Solver_* functions are thin wrappers that call generic implementations
 * from ode_interface.cpp with LorenzSystem types.
 */

#include <Rcpp.h>
#include <XAD/XAD.hpp>
#include <odelia/ode_solver.hpp>
#include <odelia/gradient.hpp>
#include <odelia/solver_interface.hpp>
#include <examples/lorenz_system.hpp>

using namespace Rcpp;
using namespace odelia;

// Define types for Lorenz system
typedef LorenzSystem<double> SystemType;
typedef LorenzSystem<xad::adj<double>::active_type> ActiveSystemType;

// Helper to get system pointer
inline Rcpp::XPtr<SystemType> get_system(SEXP xp) {
  return Rcpp::XPtr<SystemType>(xp);
}

//-------------------------------------------------------------------------
// System interface (Lorenz-specific)

// [[Rcpp::export]]
SEXP System_new(double sigma, double R, double b) {
  return Rcpp::XPtr<SystemType>(new SystemType(sigma, R, b), true);
}

// [[Rcpp::export]]
Rcpp::NumericVector System_pars(SEXP system_xp) {
  return Rcpp::wrap(get_system(system_xp)->pars());
}

// [[Rcpp::export]]
void System_set_params(SEXP system_xp, Rcpp::NumericVector params) {
  auto lor = get_system(system_xp);
  std::vector<double> tmp(params.begin(), params.end());
  lor->set_params(tmp.begin());
  lor->compute_rates();
}

// [[Rcpp::export]]
void System_set_state(SEXP system_xp, Rcpp::NumericVector y, double time) {
  auto lor = get_system(system_xp);
  if (y.size() != lor->ode_size()) {
    Rcpp::stop("State vector size mismatch");
  }
  std::vector<double> tmp(y.begin(), y.end());
  lor->set_ode_state(tmp.begin(), time);
}

// [[Rcpp::export]]
Rcpp::NumericVector System_state(SEXP system_xp) {
  auto lor = get_system(system_xp);
  std::vector<double> tmp(lor->ode_size());
  lor->ode_state(tmp.begin());
  return Rcpp::wrap(tmp);
}

// [[Rcpp::export]]
void System_set_initial_state(SEXP system_xp, Rcpp::NumericVector y, double t0 = 0.0) {
  auto lor = get_system(system_xp);
  if (y.size() != lor->ode_size()) {
    Rcpp::stop("State vector size mismatch");
  }
  std::vector<double> tmp(y.begin(), y.end());
  lor->set_initial_state(tmp.begin(), t0);
}

// [[Rcpp::export]]
void System_reset(SEXP system_xp) {
  auto lor = get_system(system_xp);
  lor->reset();
}

// [[Rcpp::export]]
Rcpp::NumericVector System_rates(SEXP system_xp) {
  auto lor = get_system(system_xp);
  std::vector<double> tmp(lor->ode_size());
  lor->ode_rates(tmp.begin());
  return Rcpp::wrap(tmp);
}

//-------------------------------------------------------------------------
// Solver interface (Lorenz-specific creation, generic operations)

// Map an R-facing method string to the solver Method enum.
static ode::Method parse_method(const std::string& method) {
  if (method == "rodas" || method == "implicit") {
    return ode::Method::rodas;
  }
  if (method == "rkck" || method == "rk45" || method == "explicit") {
    return ode::Method::rkck;
  }
  Rcpp::stop("Unknown method '" + method + "'. Use 'rkck' or 'rodas'.");
}

// Solver creation - Lorenz-specific (must know LorenzSystem type). R only ever
// holds the double solver; gradients build the active replay internally. The
// stiff `rodas` stepper is a double-path option; the active replay stays rkck
// (the Rosenbrock linear solve is not available for the AD scalar type).
// [[Rcpp::export]]
SEXP Solver_new(SEXP system_xp, SEXP control_xp, std::string method = "rkck") {
  Rcpp::XPtr<SystemType> sys(system_xp);
  Rcpp::XPtr<ode::OdeControl> ctrl(control_xp);
  SystemType sys_copy(*sys);
  auto* solver = new ode::Solver<SystemType>(sys_copy, *ctrl, parse_method(method));
  return Rcpp::XPtr<ode::Solver<SystemType>>(solver, true);
}

// All other Solver functions call the generic (double-only) templates.

// [[Rcpp::export]]
void Solver_reset(SEXP solver_xp) {
  odelia::solver::Solver_reset_impl<SystemType>(solver_xp);
}

// [[Rcpp::export]]
double Solver_time(SEXP solver_xp) {
  return odelia::solver::Solver_time_impl<SystemType>(solver_xp);
}

// [[Rcpp::export]]
Rcpp::NumericVector Solver_state(SEXP solver_xp) {
  return odelia::solver::Solver_state_impl<SystemType>(solver_xp);
}

// [[Rcpp::export]]
Rcpp::NumericVector Solver_times(SEXP solver_xp) {
  return odelia::solver::Solver_times_impl<SystemType>(solver_xp);
}

// [[Rcpp::export]]
void Solver_set_state(SEXP solver_xp, Rcpp::NumericVector y, double time) {
  odelia::solver::Solver_set_state_impl<SystemType>(solver_xp, y, time);
}

// [[Rcpp::export]]
void Solver_advance_adaptive(SEXP solver_xp, Rcpp::NumericVector times) {
  odelia::solver::Solver_advance_adaptive_impl<SystemType>(solver_xp, times);
}

// [[Rcpp::export]]
void Solver_advance_fixed(SEXP solver_xp, Rcpp::NumericVector times) {
  odelia::solver::Solver_advance_fixed_impl<SystemType>(solver_xp, times);
}

// [[Rcpp::export]]
void Solver_advance_euler(SEXP solver_xp, Rcpp::NumericVector times) {
  odelia::solver::Solver_advance_euler_impl<SystemType>(solver_xp, times);
}

// [[Rcpp::export]]
void Solver_step(SEXP solver_xp) {
  odelia::solver::Solver_step_impl<SystemType>(solver_xp);
}

// [[Rcpp::export]]
bool Solver_get_collect(SEXP solver_xp) {
  return odelia::solver::Solver_get_collect_impl<SystemType>(solver_xp);
}

// [[Rcpp::export]]
void Solver_set_collect(SEXP solver_xp, bool x) {
  odelia::solver::Solver_set_collect_impl<SystemType>(solver_xp, x);
}

// [[Rcpp::export]]
std::size_t Solver_get_history_size(SEXP solver_xp) {
  return odelia::solver::Solver_get_history_size_impl<SystemType>(solver_xp);
}

// [[Rcpp::export]]
Rcpp::DataFrame Solver_get_history_step(SEXP solver_xp, std::size_t i) {
  return odelia::solver::Solver_get_history_step_impl<SystemType>(solver_xp, i);
}

// [[Rcpp::export]]
Rcpp::List Solver_get_history(SEXP solver_xp) {
  return odelia::solver::Solver_get_history_impl<SystemType>(
    solver_xp
  );
}

// value + least-squares gradient on the double handle: the active replay is built,
// used, and destroyed inside the call. Observations are passed per call and owned
// by the functional -- the solver holds no calibration state.
// [[Rcpp::export]]
Rcpp::List Solver_value_and_gradient(SEXP solver_xp,
                                     Rcpp::NumericVector times,
                                     Rcpp::NumericMatrix observations,
                                     Rcpp::IntegerVector obs_indices,
                                     Rcpp::Nullable<Rcpp::NumericVector> ic = R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericVector> params = R_NilValue) {
  return odelia::solver::Solver_value_and_gradient_impl<SystemType, ActiveSystemType>(
    solver_xp, ic, params, times, observations, obs_indices
  );
}

//-------------------------------------------------------------------------
// An emergent-metric stand-in: the summed final state (one scalar). A pure
// reduction -- the driver replays, this reads solver.state() -- exercising the
// functional seam with something other than a loss.
struct sum_final_state {
  std::size_t codomain() const { return 1; }
  template<typename Solver>
  typename Solver::value_type operator()(Solver& solver) const {
    typename Solver::value_type total(0.0);
    for (auto const& s : solver.state()) {
      total += s;
    }
    return total;
  }
};

// DOUBLE solver handle; the active replay is built internally.
// [[Rcpp::export]]
Rcpp::List Solver_gradient_final_state(SEXP solver_xp,
                                       Rcpp::NumericVector times,
                                       Rcpp::NumericVector params) {
  std::vector<double> t(times.begin(), times.end());
  // Seed the parameters.
  ode::DifferentiationTargets targets;
  for (int i = 0; i < params.size(); ++i) {
    targets.params.push_back(i);
    targets.values.push_back(params[i]);
  }
  // `times` is the replay schedule the driver advances through.
  auto [value, gradient] = odelia::solver::gradient_on_double<SystemType, ActiveSystemType>(
    solver_xp, targets, t, sum_final_state{});
  return Rcpp::List::create(Rcpp::Named("value") = value,
                            Rcpp::Named("gradient") = Rcpp::wrap(gradient));
}

// A multi-output functional: the whole final state. Its codomain is the ODE size,
// carried so the driver reads one consistent output count.
struct final_state {
  std::size_t m;
  std::size_t codomain() const { return m; }
  template<typename Solver>
  std::vector<typename Solver::value_type> operator()(Solver& solver) const {
    return solver.state();
  }
};

// DOUBLE solver handle; the active replay is built internally.
// [[Rcpp::export]]
Rcpp::List Solver_jacobian_final_state(SEXP solver_xp,
                                       Rcpp::NumericVector times,
                                       Rcpp::NumericVector params) {
  std::vector<double> t(times.begin(), times.end());
  // Seed the parameters.
  ode::DifferentiationTargets targets;
  for (int i = 0; i < params.size(); ++i) {
    targets.params.push_back(i);
    targets.values.push_back(params[i]);
  }
  // final_state returns the whole ODE state, so its codomain is the ODE size.
  const std::size_t m =
    odelia::solver::get_solver<SystemType>(solver_xp)->get_system_ref().ode_size();
  // `times` is the replay schedule the driver advances through.
  auto [values, jacobian] = odelia::solver::jacobian_on_double<SystemType, ActiveSystemType>(
    solver_xp, targets, t, final_state{m});

  const size_t nrow = jacobian.size(), ncol = nrow ? jacobian[0].size() : 0;
  Rcpp::NumericMatrix J(nrow, ncol);
  for (size_t i = 0; i < nrow; ++i) {
    for (size_t j = 0; j < ncol; ++j) {
      J(i, j) = jacobian[i][j];
    }
  }
  return Rcpp::List::create(Rcpp::Named("values") = Rcpp::wrap(values),
                            Rcpp::Named("jacobian") = J);
}

//-------------------------------------------------------------------------
// Comparison function for deSolve

// [[Rcpp::export]]
List lorenz_rhs(double t, NumericVector state, NumericVector pars) {
  double x = state["x"], y = state["y"], z = state["z"];
  double sigma = pars["sigma"], rho = pars["rho"], beta = pars["beta"];
  
  double dx = sigma * (y - x);
  double dy = x * (rho - z) - y;
  double dz = x * y - beta * z;
  
  return List::create(NumericVector::create(dx, dy, dz));
}

//-------------------------------------------------------------------------
// Test function for parameter type conversion

// [[Rcpp::export]]
Rcpp::List test_param_types(SEXP system_xp) {
  using ad = xad::adj<double>;
  using ad_type = ad::active_type;
  using ActiveSys = LorenzSystem<ad_type>;
  
  Rcpp::XPtr<SystemType> sys(system_xp);
  
  // Get params from passive system (doubles)
  auto params = sys->pars();
  
  // Create active system using double params (should cast to ad_type)
  auto* sys_active = new ActiveSys(params[0], params[1], params[2]);
  
  // Get params back from active system
  auto params_active = sys_active->pars();
  
  // Test if gradient flows through params
  ad::tape_type tape;
  
  std::vector<ad_type> test_params(3);
  test_params[0] = params[0];
  test_params[1] = params[1]; 
  test_params[2] = params[2];
  
  tape.registerInput(test_params[0]);
  tape.newRecording();
  
  // Recreate system with registered params
  ActiveSys sys_test(test_params[0], test_params[1], test_params[2]);
  auto retrieved = sys_test.pars();
  
  ad_type result = retrieved[0] * 2.0;
  
  tape.registerOutput(result);
  xad::derivative(result) = 1.0;
  tape.computeAdjoints();
  
  double gradient = xad::derivative(test_params[0]);
  
  delete sys_active;
  
  return Rcpp::List::create(
    Rcpp::Named("passive_params") = params,
    Rcpp::Named("active_params") = params_active,
    Rcpp::Named("gradient_test") = gradient,
    Rcpp::Named("gradient_should_be") = 2.0
  );
}