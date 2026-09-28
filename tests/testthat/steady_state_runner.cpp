// Test-only steady-state / IFT-sensitivity runner, compiled on demand via
// Rcpp::sourceCpp (following vanderpol_runner.cpp from #37).
//
// The in-package Lorenz example has no fixed-point attractor at its usual
// (chaotic) parameters, so it is unsuitable for validating an equilibrium
// solver. This small autonomous system has a closed-form equilibrium and
// closed-form parameter sensitivity, so the Newton solve, the implicit-function
// -theorem sensitivity, and the eigenvalue-based stability check can all be
// checked against analytic answers and against finite differences.
//
//   y0' = a - b*y0
//   y1' = y0^2 - c*y1
//
// Equilibrium:  y0* = a/b,  y1* = a^2 / (b^2 c).
// State Jacobian df/dy = [[-b, 0], [2 y0*, -c]]  -> eigenvalues -b, -c
//   (attracting for b, c > 0). The y0^2 term makes it genuinely nonlinear, so
//   Newton takes several iterations from a poor guess.
// Parameters theta = (a, b, c); analytic sensitivity dy*/dtheta:
//   dy0*/da = 1/b            dy0*/db = -a/b^2         dy0*/dc = 0
//   dy1*/da = 2a/(b^2 c)     dy1*/db = -2 a^2/(b^3 c) dy1*/dc = -a^2/(b^2 c^2)

// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <vector>
#include <cstddef>
#include <XAD/XAD.hpp>
#include <odelia/ode_solver.hpp>
#include <odelia/ode_steady_state.hpp>

using namespace odelia;

template <typename T = double>
class DemogSystem {
public:
  using value_type = T;

  DemogSystem(T a_, T b_, T c_)
    : a(a_), b(b_), c(c_),
      y0_init(0.0), y1_init(0.0), t0(0.0),
      y0(0.0), y1(0.0), d0(0.0), d1(0.0), time(0.0) {
    reset();
  }

  size_t ode_size() const { return 2; }
  double ode_time() const { return time; }
  double ode_t0() const { return t0; }

  template <typename It>
  It set_ode_state(It it, double time_) {
    time = time_;
    y0 = *it++;
    y1 = *it++;
    compute_rates();
    return it;
  }

  void compute_rates() {
    d0 = a - b * y0;
    d1 = y0 * y0 - c * y1;
  }

  template <typename It>
  It set_initial_state(It it, double t0_ = 0.0) {
    t0 = t0_;
    y0_init = *it++;
    y1_init = *it++;
    return it;
  }

  template <typename It>
  It ode_state(It it) const { *it++ = y0; *it++ = y1; return it; }

  template <typename It>
  It ode_initial_state(It it) const {
    *it++ = y0_init; *it++ = y1_init; return it;
  }

  template <typename It>
  It ode_rates(It it) const { *it++ = d0; *it++ = d1; return it; }

  void reset() { y0 = y0_init; y1 = y1_init; time = t0; compute_rates(); }

  std::vector<double> pars() const {
    return { xad::value(a), xad::value(b), xad::value(c) };
  }

  // Pointers to the differentiable parameters, in a fixed order, so the
  // forward-mode parameter Jacobian df/dtheta can seed their tangents.
  std::vector<T*> ode_parameters() { return { &a, &b, &c }; }

  template <typename U>
  DemogSystem<U> rebind() const {
    DemogSystem<U> s(U(xad::value(a)), U(xad::value(b)), U(xad::value(c)));
    std::vector<U> init{ U(xad::value(y0_init)), U(xad::value(y1_init)) };
    s.set_initial_state(init.begin(), t0);
    std::vector<U> st{ U(xad::value(y0)), U(xad::value(y1)) };
    s.set_ode_state(st.begin(), time);
    return s;
  }

private:
  T a, b, c;
  T y0_init, y1_init;
  double t0;
  T y0, y1;
  T d0, d1;
  double time;
};

using Sys = DemogSystem<double>;

static Sys make_system(const std::vector<double>& theta) {
  return Sys(theta[0], theta[1], theta[2]);
}

// Solve to equilibrium (optionally with a RODAS warm-start through the
// transient) and report everything: the fixed point, convergence diagnostics,
// the IFT sensitivity dy*/dtheta, the eigenvalues of df/dy, and the stability
// verdict.
// [[Rcpp::export]]
Rcpp::List ss_run(std::vector<double> theta, std::vector<double> y0,
                  bool warmup) {
  Sys sys = make_system(theta);

  ode::SteadyState<Sys> ss;
  ode::SteadyState<Sys>::Options opt;
  opt.tol = 1e-12;

  ode::SteadyState<Sys>::Result res;
  if (warmup) {
    ode::OdeControl ctrl(1e-10, 1e-10, 1.0, 0.0, 1e-12, 100.0, 1e-6);
    std::vector<double> times;
    for (double t = 0.0; t <= 50.0; t += 1.0) {
      times.push_back(t);
    }
    res = ss.solve_with_warmup(sys, y0, ctrl, ode::Method::rodas, times, opt);
  } else {
    res = ss.solve(sys, y0, opt);
  }

  size_t n_params = 0;
  std::vector<double> S = ss.sensitivity(sys, n_params);
  const size_t n = res.y.size();
  Rcpp::NumericMatrix sens(n, n_params);
  for (size_t r = 0; r < n; ++r) {
    for (size_t c = 0; c < n_params; ++c) {
      sens(r, c) = S[r * n_params + c];
    }
  }

  std::vector<double> re, im;
  ss.eigenvalues(re, im);

  return Rcpp::List::create(
      Rcpp::Named("y") = res.y,
      Rcpp::Named("residual_norm") = res.residual_norm,
      Rcpp::Named("iterations") = static_cast<int>(res.iterations),
      Rcpp::Named("converged") = res.converged,
      Rcpp::Named("sensitivity") = sens,
      Rcpp::Named("eig_re") = re,
      Rcpp::Named("eig_im") = im,
      Rcpp::Named("spectral_abscissa") = ss.spectral_abscissa(),
      Rcpp::Named("stable") = ss.is_stable(),
      Rcpp::Named("time_dependence") = ss.time_dependence());
}

// Equilibrium only, for finite-difference validation of the sensitivity from R.
// [[Rcpp::export]]
std::vector<double> ss_equilibrium(std::vector<double> theta,
                                   std::vector<double> y0) {
  Sys sys = make_system(theta);
  ode::SteadyState<Sys> ss;
  ode::SteadyState<Sys>::Options opt;
  opt.tol = 1e-13;
  ode::SteadyState<Sys>::Result res = ss.solve(sys, y0, opt);
  if (!res.converged) {
    Rcpp::stop("ss_equilibrium: Newton did not converge");
  }
  return res.y;
}
