// Test-only consumer-resource System, compiled on demand via Rcpp::sourceCpp.
//
// A consumer N feeding on a resource R with a saturating (type II) uptake, the
// resource regrowing logistically (Rosenzweig & MacArthur 1963):
//   dN/dt = N (a R / (h + R) - m)
//   dR/dt = r R (1 - R / K) - q N a R / (h + R)
// A small half-saturation h makes the resource relax much faster than the
// consumer, so the system is stiff there: what RODAS is for, and the reason a
// gradient through it is wanted.

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <XAD/XAD.hpp>
#include <odelia/ode_solver.hpp>
#include <odelia/solver_interface.hpp>

using namespace odelia;

template <typename T = double>
class ConsumerResource {
public:
  using value_type = T;
  template <typename> friend class ConsumerResource;

  explicit ConsumerResource(std::vector<double> p = {1, 0.1, 0.5, 0.2, 1, 2})
    : a(p[0]), h(p[1]), q(p[2]), m(p[3]), r(p[4]), K(p[5]) {
    reset();
  }

  template <class S>
  ConsumerResource<S> rebind() const {
    ConsumerResource<S> out({xad::value(a), xad::value(h), xad::value(q),
                             xad::value(m), xad::value(r), xad::value(K)});
    out.N_init = xad::value(N_init);
    out.R_init = xad::value(R_init);
    out.t0 = t0;
    std::vector<S> y{S(xad::value(N)), S(xad::value(R))};
    out.set_ode_state(y.begin(), time);
    return out;
  }

  size_t ode_size() const { return 2; }
  double ode_time() const { return time; }

  template <class It> It set_ode_state(It it, double t) {
    N = *it++; R = *it++; time = t;
    compute_rates();
    return it;
  }
  void set_recorded_state(const std::vector<T>& y, double t) {
    set_ode_state(y.begin(), t);
  }
  template <class It> It set_initial_state(It it, double t) {
    N_init = *it++; R_init = *it++; t0 = t;
    return it;
  }
  template <class It> It set_params(It it) {
    a = *it++; h = *it++; q = *it++; m = *it++; r = *it++; K = *it++;
    return it;
  }
  void reset() {
    std::vector<T> y{N_init, R_init};
    set_ode_state(y.begin(), t0);
  }

  void compute_rates() {
    const T uptake = a * R / (h + R);
    dN = N * (uptake - m);
    dR = r * R * (T(1.0) - R / K) - q * N * uptake;
  }

  template <class It> It ode_state(It it) const { *it++ = N; *it++ = R; return it; }
  template <class It> It ode_rates(It it) const { *it++ = dN; *it++ = dR; return it; }

  std::vector<T*> ad_parameters() { return {&a, &h, &q, &m, &r, &K}; }
  template <class F> void for_each_active(F&& f) {
    f(a); f(h); f(q); f(m); f(r); f(K); f(N); f(R); f(dN); f(dR);
    f(N_init); f(R_init);
  }

private:
  T a, h, q, m, r, K;
  T N_init = T(0.5), R_init = T(1.0);
  T N = T(0.0), R = T(0.0), dN = T(0.0), dR = T(0.0);
  double t0 = 0.0, time = 0.0;
};

using CR = ConsumerResource<double>;

ode::Method method_of(const std::string& method) {
  if (method == "rodas") return ode::Method::rodas;
  if (method == "dopri") return ode::Method::dopri;
  return ode::Method::rkck;
}

ode::OdeControl control_of(double tol) {
  ode::OdeControl control;
  control.set_tol_abs(tol);
  control.set_tol_rel(tol);
  return control;
}

// A solver at parameters `p` and initial state `y0`, as R holds one.
// [[Rcpp::export]]
SEXP cr_solver(std::vector<double> p, std::vector<double> y0,
               std::string method, double tol) {
  CR sys(p);
  sys.set_initial_state(y0.begin(), 0.0);
  sys.reset();
  return Rcpp::XPtr<ode::Solver<CR>>(
      new ode::Solver<CR>(sys, control_of(tol), method_of(method)), true);
}

// [[Rcpp::export]]
Rcpp::List cr_sensitivities(SEXP solver, Rcpp::NumericVector times,
                            Rcpp::IntegerVector parameters,
                            Rcpp::IntegerVector initial) {
  return solver::Solver_sensitivities_impl<CR>(solver, times, parameters, initial);
}

// The same, two lanes at a time.
// [[Rcpp::export]]
Rcpp::List cr_sensitivities_narrow(SEXP solver, Rcpp::NumericVector times,
                                   Rcpp::IntegerVector parameters,
                                   Rcpp::IntegerVector initial) {
  return solver::Solver_sensitivities_impl<CR, 2>(solver, times, parameters,
                                                  initial);
}

// [[Rcpp::export]]
Rcpp::List cr_loss_and_gradient(SEXP solver, Rcpp::NumericVector schedule,
                                Rcpp::NumericMatrix target,
                                Rcpp::IntegerVector obs_indices,
                                Rcpp::NumericVector initial_state,
                                Rcpp::NumericVector params, std::string mode) {
  return solver::Solver_loss_and_gradient_impl<CR>(
      solver, schedule, target, obs_indices, initial_state, params, mode);
}

// The step times an adaptive run from `y0` lands on, every output time among
// them.
// [[Rcpp::export]]
std::vector<double> cr_schedule(std::vector<double> p, std::vector<double> y0,
                                std::vector<double> times, std::string method,
                                double tol) {
  CR sys(p);
  sys.set_initial_state(y0.begin(), times[0]);
  sys.reset();
  ode::Solver<CR> s(sys, control_of(tol), method_of(method));
  s.set_collect(false);
  s.advance_adaptive(times);
  return s.times();
}

// The state at every time of `schedule`, stepped through it at the sizes it
// implies: a pinned run, to take central differences of.
// [[Rcpp::export]]
Rcpp::NumericMatrix cr_pinned(std::vector<double> p, std::vector<double> y0,
                              std::vector<double> schedule, std::string method) {
  CR sys(p);
  sys.set_initial_state(y0.begin(), schedule[0]);
  sys.reset();
  ode::Solver<CR> s(sys, ode::OdeControl(), method_of(method));
  s.set_collect(false);
  std::vector<ode::schedule_row> rows{
      {schedule[0], std::numeric_limits<double>::quiet_NaN()}};
  for (size_t k = 1; k < schedule.size(); ++k) {
    rows.push_back({schedule[k], schedule[k] - schedule[k - 1]});
  }
  s.set_keep_states(true);
  s.advance_schedule(rows);
  const auto rec = s.recording();
  Rcpp::NumericMatrix out(rec.size(), 2);
  for (size_t k = 0; k < rec.size(); ++k) {
    out(k, 0) = rec[k].state[0];
    out(k, 1) = rec[k].state[1];
  }
  return out;
}
