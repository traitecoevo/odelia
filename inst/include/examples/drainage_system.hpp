#ifndef DRAINAGE_SYSTEM_HPP_
#define DRAINAGE_SYSTEM_HPP_

#include <odelia/mri.hpp>
#include <XAD/XAD.hpp>
#include <vector>
#include <cmath>

// A stiff multirate demonstrator shaped like the TF24 soil column: L fast layers
// that drain by a power law (a stiff, wet-end nonlinearity with a closed-form
// recession) and refill gently toward a slow signal, coupled to M slow modes that
// track the layer mean. The drainage is the stiff part; splitting integrates it
// exactly (analytic_flow) and leaves ROS34PW2 only the gentle refill, so the same
// model can be run with and without splitting to measure what an analytic-flow
// exposure buys as the drainage stiffness grows. The drainage rate c and the
// initial state are seedable, so gradients can flow through the exact recession.
//
// Fast:  u_l' = -c u_l^p            (drainage, stiff; recession u_l = [u_l^{1-p} + (p-1) c t]^{-1/(p-1)})
//              + r (g - u_l)        (gentle refill toward the slow signal g)
// Slow:  x_j' = w_j (mean(u) - x_j)
// g = mean(x).
template <typename T = double>
class DrainageSystem {
public:
  using value_type = T;
  using vec = std::vector<T>;

  DrainageSystem(double drain_c, int n_fast_, int n_slow_)
    : c(drain_c), p(16.0), r(0.3), n_fast(n_fast_), n_slow(n_slow_),
      omega(n_slow_), u_init(n_fast_), x_init(n_slow_), u(n_fast_), x(n_slow_),
      t0(0.0), time(0.0) {
    for (int j = 0; j < n_slow; ++j)
      omega[j] = 0.02 + 0.18 * (n_slow > 1 ? double(j) / (n_slow - 1) : 0.0);
    for (int l = 0; l < n_fast; ++l) u_init[l] = T(0.9);   // wet: drainage active
    for (int j = 0; j < n_slow; ++j) x_init[j] = T(0.5 + 0.2 * std::cos(3.14159265358979 * j / n_slow));
    reset();
  }

  size_t fast_size() const { return (size_t)n_fast; }
  size_t slow_size() const { return (size_t)n_slow; }
  size_t coupling_size() const { return 1; }
  size_t ode_size() const { return fast_size() + slow_size(); }
  double ode_time() const { return time; }
  double ode_t0() const { return t0; }

  void aggregate(const vec& xs, vec& g) const {
    T s = T(0.0);
    for (const auto& xj : xs) s += xj;
    g[0] = s / double(n_slow);
  }

  void fast_rates(const vec& us, const vec& g, vec& du) const {
    for (int l = 0; l < n_fast; ++l)
      du[l] = -drainage(us[l]) + r * (g[0] - us[l]);
  }

  void slow_rates(const vec& xs, const vec& us, vec& dx) const {
    T m = T(0.0);
    for (const auto& ul : us) m += ul;
    m /= double(n_fast);
    for (int j = 0; j < n_slow; ++j) dx[j] = omega[j] * (m - xs[j]);
  }

  // Split hooks: the exact drainage recession (at the state scalar, so it tapes),
  // and the gentle remainder (templated on the argument scalar so the value-
  // Jacobian can evaluate it at double while the stages run active).
  void analytic_flow(vec& us, double dt) const {
    using std::pow; using xad::pow;
    for (auto& ul : us) {
      if (xad::value(ul) <= 0.0) continue;
      const T base = pow(ul, 1.0 - p) + (p - 1.0) * c * dt;
      ul = pow(base, -1.0 / (p - 1.0));
    }
  }
  template <class U>
  void residual_rhs(const std::vector<U>& us, const std::vector<U>& g, std::vector<U>& du) const {
    for (int l = 0; l < n_fast; ++l) du[l] = r * (g[0] - us[l]);
  }

  // State laid out [slow x | fast u] (slow-first, the multirate layout).
  template <typename Iterator>
  Iterator set_ode_state(Iterator it, double time_) {
    time = time_;
    for (auto& xj : x) xj = *it++;
    for (auto& ul : u) ul = *it++;
    return it;
  }
  template <typename Iterator>
  Iterator ode_state(Iterator it) const {
    for (const auto& xj : x) *it++ = xj;
    for (const auto& ul : u) *it++ = ul;
    return it;
  }
  template <typename Iterator>
  Iterator ode_rates(Iterator it) const {
    vec g(1), du(n_fast), dx(n_slow);
    aggregate(x, g);
    fast_rates(u, g, du);
    slow_rates(x, u, dx);
    for (const auto& d : dx) *it++ = d;
    for (const auto& d : du) *it++ = d;
    return it;
  }

  template <typename Iterator>
  Iterator set_params(Iterator it) { c = *it++; return it; }
  template <typename Tape, typename Iterator>
  std::vector<T*> set_params(Tape& tape, Iterator it) {
    c = *it++;
    tape.registerInput(c);
    return {&c};
  }

  // Initial state laid out [slow x | fast u], matching ode_state.
  template <typename Iterator>
  Iterator set_initial_state(Iterator it, double t0_ = 0.0) {
    t0 = t0_;
    for (auto& xj : x_init) xj = *it++;
    for (auto& ul : u_init) ul = *it++;
    return it;
  }
  template <typename Tape, typename Iterator>
  std::vector<T*> set_initial_state(Tape& tape, Iterator it, double t0_ = 0.0) {
    t0 = t0_;
    std::vector<T*> refs;
    for (auto& xj : x_init) { xj = *it++; tape.registerInput(xj); refs.push_back(&xj); }
    for (auto& ul : u_init) { ul = *it++; tape.registerInput(ul); refs.push_back(&ul); }
    return refs;
  }

  void reset() {
    for (int l = 0; l < n_fast; ++l) u[l] = u_init[l];
    for (int j = 0; j < n_slow; ++j) x[j] = x_init[j];
    time = t0;
  }

  std::vector<double> pars() const { return {xad::value(c), double(n_fast), double(n_slow)}; }

private:
  T drainage(const T& ul) const { using std::pow; using xad::pow; return c * pow(ul, p); }

  T c;
  double p, r;
  int n_fast, n_slow;
  std::vector<double> omega;
  vec u_init, x_init, u, x;
  double t0, time;
};

#endif
