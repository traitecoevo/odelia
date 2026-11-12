// The times a sweep or a replay hands a System are the times the run handed it.
//
// A System that reads a record of another run by time -- a field recorded per
// stage, looked up at the stage's time -- finds nothing if a stage time is an ulp
// from the run's. Each is `time + c * h` on the same doubles, so they can differ
// only where a compiler rounds two copies differently, fusing the multiply-add in
// one and not the other. So every copy is held to one rounding, the product and
// then the sum. This is built with optimisation and the compiler's own
// contraction default, which is where a copy could be fused.

#include <odelia/ode_solver.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

int failures = 0;

void check(bool ok, const std::string& what) {
  std::printf("  %s %s\n", ok ? "ok  " : "FAIL", what.c_str());
  if (!ok) {
    ++failures;
  }
}

// Every time a System is set to, at either scalar, in call order.
std::vector<double> handed;

// A damped spring, x'' = -k x - c x', whose period is short against the run, so
// the run takes thousands of steps of varying size.
template <typename T = double>
struct Spring {
  using value_type = T;
  template <typename> friend struct Spring;
  T k, x, v, dx, dv, x_init, v_init;
  double time = 0.0, t0 = 0.0;
  explicit Spring(T k_ = T(400.0))
      : k(k_), x(1.0), v(0.0), dx(0.0), dv(0.0), x_init(1.0), v_init(0.0) {
    compute_rates();
  }
  template <class S> Spring<S> rebind() const {
    Spring<S> s(S(xad::value(k)));
    s.x_init = S(xad::value(x_init));
    s.v_init = S(xad::value(v_init));
    s.x = S(xad::value(x));
    s.v = S(xad::value(v));
    s.time = time;
    s.t0 = t0;
    s.compute_rates();
    return s;
  }
  size_t ode_size() const { return 2; }
  double ode_time() const { return time; }
  template <class It> It set_ode_state(It it, double t) {
    x = *it++;
    v = *it++;
    time = t;
    handed.push_back(t);
    compute_rates();
    return it;
  }
  void set_recorded_state(const std::vector<T>& s, double t) {
    set_ode_state(s.begin(), t);
  }
  void compute_rates() {
    dx = v;
    dv = -k * x - T(0.5) * v;
  }
  template <class It> It ode_state(It it) const {
    *it++ = x;
    *it++ = v;
    return it;
  }
  template <class It> It ode_rates(It it) const {
    *it++ = dx;
    *it++ = dv;
    return it;
  }
  template <class It> It ode_initial_state(It it) const {
    *it++ = x_init;
    *it++ = v_init;
    return it;
  }
  template <class It> It set_initial_state(It it, double t0_ = 0.0) {
    t0 = t0_;
    x_init = *it++;
    v_init = *it++;
    return it;
  }
  void reset() {
    x = x_init;
    v = v_init;
    time = t0;
    compute_rates();
  }
  std::vector<T*> ad_parameters() { return {&k}; }
  template <class F> void for_each_active(F&& f) {
    f(k); f(x); f(v); f(dx); f(dv); f(x_init); f(v_init);
  }
};

// The interior stage times of a step of `h` from `start`, as the steppers round
// them: the product on its own, then the sum.
std::vector<double> stage_times(double start, double h,
                                const std::vector<double>& fractions) {
  std::vector<double> out;
  for (const double c : fractions) {
    const volatile double offset = c * h;
    out.push_back(start + offset);
  }
  return out;
}

// How many of `times` are absent from the sorted `among`, and the first of them.
size_t count_absent(const std::vector<double>& times,
                    const std::vector<double>& among, double& first) {
  size_t absent = 0;
  for (const double t : times) {
    if (!std::binary_search(among.begin(), among.end(), t) && absent++ == 0) {
      first = t;
    }
  }
  return absent;
}

void report(const char* name, const char* what, size_t absent, double first) {
  std::printf("       %s: %zu %s", name, absent, what);
  if (absent > 0) {
    std::printf(", first %a", first);
  }
  std::printf("\n");
}

void test_sweep_hands_the_run_times(odelia::ode::Method m, const char* name,
                                    const std::vector<double>& fractions) {
  Spring<double> sys;
  odelia::ode::OdeControl control(1e-10, 1e-10, 1.0, 0.0, 1e-12, 1.0, 1e-4);
  odelia::ode::Solver<Spring<double>> s(sys, control, m);
  s.set_keep_states(true);
  handed.clear();
  s.reset();
  s.advance_adaptive(std::vector<double>{0.0, 20.0});
  // The run stands at each row's time, its start among them, whether or not a
  // rate evaluation handed the System that time.
  std::vector<double> run = handed;
  const auto rec = s.recording();
  for (const auto& row : rec) {
    run.push_back(row.time);
  }
  std::sort(run.begin(), run.end());
  const size_t steps = rec.size() - 1;
  check(steps > 1000, std::string(name) + ": the run takes many steps");

  // Every accepted step's interior stages, rounded as the steppers round them.
  std::vector<double> expected;
  for (size_t k = 1; k < rec.size(); ++k) {
    const auto t = stage_times(rec[k - 1].time, rec[k].step_size, fractions);
    expected.insert(expected.end(), t.begin(), t.end());
  }

  handed.clear();
  Spring<double> again;
  odelia::ode::Solver<Spring<double>> r(again, control, m);
  r.reset();
  r.replay(rec);
  std::vector<double> replayed = handed;
  std::sort(replayed.begin(), replayed.end());

  handed.clear();
  odelia::ode::adjoint_batch lambda = odelia::ode::adjoint_batch::one({1.0, 0.0});
  odelia::ode::adjoint_batch dp(1, 1);
  s.sweep_adjoint(lambda, dp);
  std::vector<double> swept = handed;
  std::sort(swept.begin(), swept.end());
  check(swept.size() >= fractions.size() * steps,
        std::string(name) + ": the sweep hands the System every stage of every step");
  check(replayed.size() >= fractions.size() * steps,
        std::string(name) + ": and so does a replay");

  double first = NAN;
  size_t absent = count_absent(expected, run, first);
  report(name, "stage times absent from the run", absent, first);
  check(absent == 0, std::string(name) + ": the run rounds each stage time as the product, then the sum");
  absent = count_absent(expected, replayed, first);
  report(name, "absent from its replay", absent, first);
  check(absent == 0, std::string(name) + ": so does a replay of it");
  absent = count_absent(expected, swept, first);
  report(name, "absent from its sweep", absent, first);
  check(absent == 0, std::string(name) + ": and so does its sweep");
  absent = count_absent(handed, run, first);
  report(name, "swept times the run never handed", absent, first);
  check(absent == 0,
        std::string(name) + ": every time the sweep hands the System, the run handed it");
}

} // namespace

int main() {
  std::printf("stage times\n");
  // The tableaux's interior stage fractions: Cash-Karp's and Dormand-Prince's.
  test_sweep_hands_the_run_times(odelia::ode::Method::rkck, "rkck",
                                 {0.2, 0.3, 0.6, 1.0, 0.875});
  test_sweep_hands_the_run_times(odelia::ode::Method::dopri, "dopri",
                                 {0.2, 0.3, 0.8, 8.0 / 9.0, 1.0});
  if (failures > 0) {
    std::printf("%d check(s) failed\n", failures);
    return 1;
  }
  std::printf("all stage-time checks passed\n");
  return 0;
}
