// A consumer-resource model, solved and differentiated from one file.
//
//   dN/dt = N (a R / (h + R) - m)
//   dR/dt = r R (1 - R / K) - q N a R / (h + R)
//
// A consumer N takes up a resource R with a saturating (type II) response,
// a R / (h + R), and is lost at a rate m, harvesting or mortality; the resource
// regrows logistically to a capacity K, and q is the resource consumed per
// consumer produced (Rosenzweig & MacArthur 1963, Am. Nat. 97:209). A small
// half-saturation h makes the resource relax much faster than the consumer,
// which is where the implicit stepper earns its place.
//
// Build it from R with odelia loaded first, so the compiled tape runtime is
// resolved from odelia's own library (ARCHITECTURE.md):
//
//   library(odelia)
//   Sys.setenv(PKG_CPPFLAGS = paste0(
//     "-I", shQuote(system.file("include", package = "odelia")),
//     " -DXAD_NO_THREADLOCAL -DXAD_USE_STRONG_INLINE"),
//     PKG_LIBS = shQuote(odelia_load_dll()))
//   Rcpp::sourceCpp(system.file("examples/consumer_resource/consumer_resource.cpp",
//                               package = "odelia"))
//
// Every scalar function the rates call must be found by unqualified lookup, so
// that XAD's overloads are chosen at an AD scalar: write `using std::exp;` and
// then `exp(x)`, never `std::exp(x)`. A quantity that must stay positive, or a
// state whose components differ by orders of magnitude, is better carried as
// its logarithm: the stepper then cannot cross zero, and the tolerances mean
// the same thing for every component.

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <string>
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
  // The rates compile at a tangent of a tangent, which a forward pass over a
  // RODAS run takes the Jacobian at.
  static constexpr bool nests_tangents = true;

  // a, h, q, m, r, K
  explicit ConsumerResource(std::vector<double> p = {1, 0.5, 0.5, 0.4, 1, 1})
    : a(p[0]), h(p[1]), q(p[2]), m(p[3]), r(p[4]), K(p[5]) {
    reset();
  }

  // --- For solving ------------------------------------------------------------
  size_t ode_size() const { return 2; }
  double ode_time() const { return time; }

  template <class It> It set_ode_state(It it, double t) {
    N = *it++; R = *it++; time = t;
    compute_rates();
    return it;
  }
  template <class It> It ode_state(It it) const { *it++ = N; *it++ = R; return it; }
  template <class It> It ode_rates(It it) const { *it++ = dN; *it++ = dR; return it; }

  void compute_rates() {
    const T uptake = a * R / (h + R);
    dN = N * (uptake - m);
    dR = r * R * (T(1.0) - R / K) - q * N * uptake;
  }

  template <class It> It set_initial_state(It it, double t) {
    N_init = *it++; R_init = *it++; t0 = t;
    return it;
  }
  void reset() {
    std::vector<T> y{N_init, R_init};
    set_ode_state(y.begin(), t0);
  }

  // --- For a derivative (Sweepable and TangentSweepable, ode_interface.hpp) ---
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
  void set_recorded_state(const std::vector<T>& y, double t) {
    set_ode_state(y.begin(), t);
  }
  std::vector<T*> ad_parameters() { return {&a, &h, &q, &m, &r, &K}; }
  template <class F> void for_each_active(F&& f) {
    f(a); f(h); f(q); f(m); f(r); f(K);
    f(N); f(R); f(dN); f(dR); f(N_init); f(R_init);
  }

private:
  T a, h, q, m, r, K;
  T N_init = T(0.5), R_init = T(0.5);
  T N = T(0.0), R = T(0.0), dN = T(0.0), dR = T(0.0);
  double t0 = 0.0, time = 0.0;
};

using CR = ConsumerResource<double>;

ode::Method method_of(const std::string& method) {
  if (method == "rodas") return ode::Method::rodas;
  if (method == "rkck") return ode::Method::rkck;
  return ode::Method::dopri;
}

ode::OdeControl control_of(double tol) {
  ode::OdeControl control;
  control.set_tol_abs(tol);
  control.set_tol_rel(tol);
  return control;
}

// A solver at parameters p and initial state y0, from time times[0], for R to
// hold: $sensitivities() and the like are instantiated on it below.
// [[Rcpp::export]]
SEXP cr_solver(std::vector<double> p, std::vector<double> y0, double t0 = 0.0,
               std::string method = "dopri", double tol = 1e-8) {
  CR system(p);
  system.set_initial_state(y0.begin(), t0);
  system.reset();
  return Rcpp::XPtr<ode::Solver<CR>>(
      new ode::Solver<CR>(system, control_of(tol), method_of(method)), true);
}

// The state at each of `times`, the first of which is where the solver stands.
// [[Rcpp::export]]
Rcpp::NumericMatrix cr_solve(SEXP solver, std::vector<double> times) {
  ode::Solver<CR> s(*solver::get_solver<CR>(solver));
  s.set_collect(false);
  Rcpp::NumericMatrix out(times.size(), 3);
  for (size_t k = 0; k < times.size(); ++k) {
    if (k > 0) s.advance_adaptive(std::vector<double>{times[k - 1], times[k]});
    out(k, 0) = times[k];
    out(k, 1) = s.state()[0];
    out(k, 2) = s.state()[1];
  }
  Rcpp::colnames(out) = Rcpp::CharacterVector::create("time", "N", "R");
  return out;
}

// y and dy/d(parameters, initial state) at each of `times`, by forward mode,
// under any stepper. The inputs are named by their place in ad_parameters()
// (a, h, q, m, r, K) and in the state (N, R).
// [[Rcpp::export]]
Rcpp::List cr_sensitivities(SEXP solver, Rcpp::NumericVector times,
                            Rcpp::IntegerVector parameters,
                            Rcpp::IntegerVector initial) {
  return solver::Solver_sensitivities_impl<CR>(solver, times, parameters, initial);
}
