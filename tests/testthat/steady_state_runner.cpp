// Test-only steady-state / IFT-sensitivity runner, compiled on demand via
// Rcpp::sourceCpp (following vanderpol_runner.cpp from #37).
//
// The in-package Lorenz example has no fixed-point attractor at its usual
// (chaotic) parameters, so it is unsuitable for validating an equilibrium
// solver. Three small autonomous systems stand in:
//
// DemogSystem, with a closed-form equilibrium and closed-form parameter
// sensitivity, so the Newton solve, the implicit-function-theorem sensitivity,
// and the eigenvalue-based stability check can all be checked against analytic
// answers and against finite differences.
//
//   y0' = a - b*y0
//   y1' = y0^2 - c*y1
//
// Equilibrium:  y0* = a/b,  y1* = a^2 / (b^2 c).
// State Jacobian df/dy = [[-b, 0], [2 y0*, -c]]  -> eigenvalues -b, -c
//   (attracting for b, c > 0; repelling in y0 for b < 0, which Newton reaches
//   just as readily). The y0^2 term makes it genuinely nonlinear, so Newton
//   takes several iterations from a poor guess.
// Parameters theta = (a, b, c); analytic sensitivity dy*/dtheta:
//   dy0*/da = 1/b            dy0*/db = -a/b^2         dy0*/dc = 0
//   dy1*/da = 2a/(b^2 c)     dy1*/db = -2 a^2/(b^3 c) dy1*/dc = -a^2/(b^2 c^2)
//
// LogisticSystem, with two roots: the trivial one at y = 0, repelling when
// r > m, and the interior attractor. Newton from a small guess lands on the
// trivial root; this is what solve_with_warmup() has to get past.
//
//   y0' = r*y0*(1 - y0/K) - m*y0
//   y1' = y0^2 - c*y1
//
// Interior equilibrium:  y0* = K (1 - m/r),  y1* = y0*^2 / c.
// df/dy there = [[m - r, 0], [2 y0*, -c]] -> eigenvalues m - r, -c.
// Parameters theta = (r, K, m, c):
//   dy0*/dr = K m / r^2   dy0*/dK = 1 - m/r   dy0*/dm = -K/r   dy0*/dc = 0
//   dy1*/d(r,K,m) = (2 y0*/c) dy0*/d(r,K,m),   dy1*/dc = -y0*^2 / c^2
//
// CachedSystem, which breaks the ad_parameters() contract on purpose: its rate
// reads k = a*b, computed once in the constructor, so seeding a or b in place
// reaches nothing. y' = k - y; y* = a*b; dy*/da = b, dy*/db = a, which the AD
// sensitivity reports as 0 and 0. check_parameters() must say so.

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
  // forward-mode parameter Jacobian df/dtheta can seed their tangents (and a
  // reverse-mode sweep can accumulate adjoints for them: the same hook, #59).
  std::vector<T*> ad_parameters() { return { &a, &b, &c }; }

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

template <typename T = double>
class LogisticSystem {
public:
  using value_type = T;

  LogisticSystem(T r_, T K_, T m_, T c_)
    : r(r_), K(K_), m(m_), c(c_),
      y0_init(0.0), y1_init(0.0), t0(0.0),
      y0(0.0), y1(0.0), d0(0.0), d1(0.0), time(0.0) {
    reset();
  }

  size_t ode_size() const { return 2; }
  double ode_time() const { return time; }
  double ode_t0() const { return t0; }
  bool ode_autonomous() const { return true; }

  template <typename It>
  It set_ode_state(It it, double time_) {
    time = time_;
    y0 = *it++;
    y1 = *it++;
    compute_rates();
    return it;
  }

  void compute_rates() {
    d0 = r * y0 * (1.0 - y0 / K) - m * y0;
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

  std::vector<T*> ad_parameters() { return { &r, &K, &m, &c }; }

  template <typename U>
  LogisticSystem<U> rebind() const {
    LogisticSystem<U> s(U(xad::value(r)), U(xad::value(K)), U(xad::value(m)),
                        U(xad::value(c)));
    std::vector<U> init{ U(xad::value(y0_init)), U(xad::value(y1_init)) };
    s.set_initial_state(init.begin(), t0);
    std::vector<U> st{ U(xad::value(y0)), U(xad::value(y1)) };
    s.set_ode_state(st.begin(), time);
    return s;
  }

private:
  T r, K, m, c;
  T y0_init, y1_init;
  double t0;
  T y0, y1;
  T d0, d1;
  double time;
};

template <typename T = double>
class CachedSystem {
public:
  using value_type = T;

  CachedSystem(T a_, T b_) : a(a_), b(b_), k(a_ * b_), y(0.0), d(0.0) {
    d = k - y;
  }

  size_t ode_size() const { return 1; }
  double ode_time() const { return 0.0; }
  bool ode_autonomous() const { return true; }

  template <typename It>
  It set_ode_state(It it, double) { y = *it++; d = k - y; return it; }

  template <typename It>
  It ode_state(It it) const { *it++ = y; return it; }

  template <typename It>
  It ode_rates(It it) const { *it++ = d; return it; }

  std::vector<T*> ad_parameters() { return { &a, &b }; }

  template <typename U>
  CachedSystem<U> rebind() const {
    CachedSystem<U> s(U(xad::value(a)), U(xad::value(b)));
    std::vector<U> st{ U(xad::value(y)) };
    s.set_ode_state(st.begin(), 0.0);
    return s;
  }

private:
  T a, b, k;
  T y, d;
};

// Run a SteadyState on `sys` and report everything: the fixed point,
// convergence diagnostics, the IFT sensitivity dy*/dtheta, the eigenvalues of
// df/dy, the stability verdict, and the ad_parameters() self-check.
template <typename Sys>
static Rcpp::List report(Sys& sys, const std::vector<double>& y0, bool warmup,
                         int max_iter) {
  ode::SteadyState<Sys> ss;
  typename ode::SteadyState<Sys>::Options opt;
  opt.tol = 1e-12;
  opt.max_iter = static_cast<size_t>(max_iter);

  typename ode::SteadyState<Sys>::Result res;
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
      Rcpp::Named("warmed") = res.warmed,
      Rcpp::Named("sensitivity") = sens,
      Rcpp::Named("eig_re") = re,
      Rcpp::Named("eig_im") = im,
      Rcpp::Named("spectral_abscissa") = ss.spectral_abscissa(),
      Rcpp::Named("stable") = ss.is_stable(),
      Rcpp::Named("time_dependence") = ss.time_dependence(),
      Rcpp::Named("param_check") = ss.check_parameters(sys));
}

// DemogSystem, theta = (a, b, c).
// [[Rcpp::export]]
Rcpp::List ss_run(std::vector<double> theta, std::vector<double> y0,
                  bool warmup, int max_iter = 100) {
  DemogSystem<double> sys(theta[0], theta[1], theta[2]);
  return report(sys, y0, warmup, max_iter);
}

// LogisticSystem, theta = (r, K, m, c).
// [[Rcpp::export]]
Rcpp::List ss_run_logistic(std::vector<double> theta, std::vector<double> y0,
                           bool warmup, int max_iter = 100) {
  LogisticSystem<double> sys(theta[0], theta[1], theta[2], theta[3]);
  return report(sys, y0, warmup, max_iter);
}

// CachedSystem, theta = (a, b).
// [[Rcpp::export]]
Rcpp::List ss_run_cached(std::vector<double> theta, std::vector<double> y0) {
  CachedSystem<double> sys(theta[0], theta[1]);
  return report(sys, y0, false, 100);
}

// Equilibrium only, for finite-difference validation of the sensitivity from R.
// [[Rcpp::export]]
std::vector<double> ss_equilibrium(std::vector<double> theta,
                                   std::vector<double> y0) {
  DemogSystem<double> sys(theta[0], theta[1], theta[2]);
  ode::SteadyState<DemogSystem<double>> ss;
  ode::SteadyState<DemogSystem<double>>::Options opt;
  opt.tol = 1e-13;
  auto res = ss.solve(sys, y0, opt);
  if (!res.converged) {
    Rcpp::stop("ss_equilibrium: Newton did not converge");
  }
  return res.y;
}

// Eigenvalues of a square matrix, straight from ode_linalg.hpp, so the
// Hessenberg reduction and the QR sweep can be checked on matrices larger than
// the 2x2 Jacobians above (for which neither runs).
// [[Rcpp::export]]
Rcpp::List ss_eigenvalues(Rcpp::NumericMatrix A) {
  const size_t n = static_cast<size_t>(A.nrow());
  if (static_cast<size_t>(A.ncol()) != n) {
    Rcpp::stop("ss_eigenvalues: matrix must be square");
  }
  std::vector<double> a(n * n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      a[i * n + j] = A(i, j);
    }
  }
  std::vector<double> re, im;
  ode::linalg::eigenvalues(a, n, re, im);
  return Rcpp::List::create(Rcpp::Named("re") = re, Rcpp::Named("im") = im);
}
