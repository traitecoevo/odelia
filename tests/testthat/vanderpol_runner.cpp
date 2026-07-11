// Test-only Van der Pol runner, compiled on demand via Rcpp::sourceCpp.
//
// Van der Pol is a classic stiff test problem; the Lorenz example shipped with
// the package is non-stiff, so this small standalone system exercises the
// implicit RODAS stepper on a genuinely stiff problem and lets the test compare
// against deSolve. Only the passive (double) path is used here.

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <cstddef>
#include <XAD/XAD.hpp>
#include <odelia/ode_solver.hpp>

using namespace odelia;

// Van der Pol in the stiff (singular-perturbation) form:
//   y0' = y1
//   y1' = ((1 - y0^2) * y1 - y0) / eps
// small eps => stiff. Single parameter eps, templated scalar type with a
// rebind_from() lift so the implicit stepper can differentiate the RHS.
template <typename T = double>
class VdpSystem {
public:
  using value_type = T;

  VdpSystem(T eps_)
    : eps(eps_), y0_init(2.0), y1_init(0.0), t0(0.0),
      y0(2.0), y1(0.0), d0(0.0), d1(0.0), time(0.0) {
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
    d0 = y1;
    d1 = ((1.0 - y0 * y0) * y1 - y0) / eps;
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

  std::vector<double> pars() const { return { xad::value(eps) }; }

  template <class S2> using rebind = VdpSystem<S2>;

  template <typename U>
  rebind<U> rebind_from() const {
    VdpSystem<U> s(U(xad::value(eps)));
    std::vector<U> init{ U(xad::value(y0_init)), U(xad::value(y1_init)) };
    s.set_initial_state(init.begin(), t0);
    std::vector<U> st{ U(xad::value(y0)), U(xad::value(y1)) };
    s.set_ode_state(st.begin(), time);
    return s;
  }

private:
  T eps;
  T y0_init, y1_init;
  double t0;
  T y0, y1;
  T d0, d1;
  double time;
};

// Integrate Van der Pol adaptively with the requested method, collecting the
// state at each supplied output time. Returns the trajectory and the total
// number of accepted internal steps (a proxy for solver effort).
// [[Rcpp::export]]
Rcpp::List vdp_run(double eps, std::vector<double> out_times,
                   std::vector<double> y0, double tol_rel, double tol_abs,
                   std::string method) {
  using Sys = VdpSystem<double>;
  const ode::Method m =
      (method == "rodas") ? ode::Method::rodas : ode::Method::rkck;

  Sys sys(eps);
  sys.set_initial_state(y0.begin(), out_times.front());

  ode::OdeControl ctrl(tol_abs, tol_rel, 1.0, 0.0, 1e-12, 100.0, 1e-6);
  ode::Solver<Sys> solver(sys, ctrl, m);
  solver.set_collect(true);
  solver.advance_adaptive(out_times);

  auto hist = solver.get_history();
  Rcpp::NumericMatrix out(hist.size(), 2);
  for (size_t i = 0; i < hist.size(); ++i) {
    std::vector<double> s(2);
    hist[i].ode_state(s.begin());
    out(i, 0) = s[0];
    out(i, 1) = s[1];
  }

  return Rcpp::List::create(
      Rcpp::Named("state") = out,
      Rcpp::Named("n_steps") = static_cast<int>(solver.times().size()));
}
