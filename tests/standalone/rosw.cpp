// The Rosenbrock-W stepper, ROS34PW2 (ode_step_rosw.hpp): its order with an
// exact Jacobian, a Jacobian kept for a whole run and a wrong one; how many
// Jacobians it forms against RODAS on two stiff problems; when it forms one;
// and what it refuses.

#include <odelia/ode_solver.hpp>
#include <odelia/ode_callback_system.hpp>
#include <examples/lorenz_system.hpp>

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

namespace ode = odelia::ode;
using ode::CallbackSystem;
using ode::Method;
using State = std::vector<double>;

int failures = 0;

void check(bool ok, const std::string& what) {
  std::printf("  %s %s\n", ok ? "ok  " : "FAIL", what.c_str());
  if (!ok) {
    ++failures;
  }
}

template <class F>
bool refused(F&& f) {
  try {
    f();
  } catch (const std::runtime_error&) {
    return true;
  }
  return false;
}

// A driven, damped pendulum: nonlinear, with a Jacobian that moves with the
// state, and not autonomous, so df/dt enters every stage.
void pendulum(double t, const State& y, State& f) {
  f[0] = y[1];
  f[1] = -std::sin(y[0]) - 0.1 * y[1] + 0.5 * std::cos(t);
}
void pendulum_jac(double, const State& y, const State&, State& J) {
  J = {0.0, 1.0, -std::cos(y[0]), -0.1};
}
// The same, off by a fixed matrix as large as the Jacobian itself.
void pendulum_jac_wrong(double t, const State& y, const State& f, State& J) {
  pendulum_jac(t, y, f, J);
  J[0] += 0.5;
  J[1] -= 0.3;
  J[2] += 0.2;
  J[3] += 0.7;
}

// The state at t = 2 after n equal steps.
State pendulum_fixed(Method m, CallbackSystem::jac_type jac, std::size_t age,
                     int n) {
  CallbackSystem sys(pendulum, {1.0, 0.0}, 0.0, jac);
  ode::OdeControl control;
  control.set_jacobian_max_age(age);
  ode::Solver<CallbackSystem> s(sys, control, m);
  s.set_collect(false);
  std::vector<double> times(n + 1);
  for (int i = 0; i <= n; ++i) {
    times[i] = 2.0 * i / n;
  }
  s.advance_fixed(times);
  return s.state();
}

// The observed order between n and 2n steps, the smallest over 20..160.
double observed_order(Method m, CallbackSystem::jac_type jac, std::size_t age,
                      const State& exact) {
  double worst = 99.0, previous = 0.0;
  for (int n : {20, 40, 80, 160}) {
    const State y = pendulum_fixed(m, jac, age, n);
    const double err = std::hypot(y[0] - exact[0], y[1] - exact[1]);
    if (previous > 0.0) {
      worst = std::min(worst, std::log2(previous / err));
    }
    previous = err;
  }
  return worst;
}

// The order conditions of a W-method hold for any matrix in place of J, so
// the order is 3 with J formed at every step, with J formed once at the start
// and kept for the run, and with a J that is wrong at every step. RODAS, whose
// order conditions use the exact J, is the control: with the wrong one it
// falls to order 1, which is what shows the test can see a lost order.
void test_order() {
  std::printf("order with an exact, a kept and a wrong Jacobian\n");
  State exact;
  {
    CallbackSystem sys(pendulum, {1.0, 0.0}, 0.0, pendulum_jac);
    ode::OdeControl tight(1e-14, 1e-14, 1.0, 0.0, 1e-14, 0.01, 1e-6);
    ode::Solver<CallbackSystem> s(sys, tight, Method::dopri);
    s.set_collect(false);
    s.advance_adaptive({0.0, 2.0});
    exact = s.state();
  }
  const double every = observed_order(Method::ros34pw2, pendulum_jac, 1, exact);
  const double kept = observed_order(Method::ros34pw2, pendulum_jac, 0, exact);
  const double wrong =
    observed_order(Method::ros34pw2, pendulum_jac_wrong, 1, exact);
  const double rodas_wrong =
    observed_order(Method::rodas, pendulum_jac_wrong, 1, exact);
  std::printf("       (orders: every step %.2f, kept %.2f, wrong %.2f; "
              "RODAS wrong %.2f)\n", every, kept, wrong, rodas_wrong);
  check(every > 2.9, "order 3 with J formed at every step");
  check(kept > 2.9, "order 3 with J formed once and kept for the run");
  check(wrong > 2.9, "order 3 with a wrong J");
  check(rodas_wrong < 1.5, "while RODAS with the wrong J loses its order");
}

// Robertson's kinetics, and Van der Pol at eps = 1e-3.
void robertson(double, const State& y, State& f) {
  f[0] = -0.04 * y[0] + 1e4 * y[1] * y[2];
  f[2] = 3e7 * y[1] * y[1];
  f[1] = -f[0] - f[2];
}
void robertson_jac(double, const State& y, const State&, State& J) {
  J = {-0.04, 1e4 * y[2], 1e4 * y[1],
       0.04, -1e4 * y[2] - 6e7 * y[1], -1e4 * y[1],
       0.0, 6e7 * y[1], 0.0};
}
void vanderpol(double, const State& y, State& f) {
  f[0] = y[1];
  f[1] = ((1.0 - y[0] * y[0]) * y[1] - y[0]) / 1e-3;
}
void vanderpol_jac(double, const State& y, const State&, State& J) {
  J = {0.0, 1.0, (-2.0 * y[0] * y[1] - 1.0) / 1e-3, (1.0 - y[0] * y[0]) / 1e-3};
}

struct problem {
  const char* name;
  CallbackSystem::rhs_type rhs;
  CallbackSystem::jac_type jac;
  State y0;
  double t_end;
  double abs_per_rel; // Robertson's middle species peaks near 4e-5
};

struct outcome {
  State y;
  std::size_t accepted, rejected, jac, rhs;
};

outcome solve(const problem& p, Method m, double tol, std::size_t age) {
  CallbackSystem sys(p.rhs, p.y0, 0.0, p.jac, CallbackSystem::valid_type(),
                     true);
  ode::OdeControl control(tol * p.abs_per_rel, tol, 1.0, 0.0, 1e-14, p.t_end,
                          1e-6);
  control.set_jacobian_max_age(age);
  ode::Solver<CallbackSystem> s(sys, control, m);
  s.set_collect(false);
  s.advance_adaptive({0.0, p.t_end});
  const CallbackSystem& after = s.get_system_ref();
  return {s.state(), s.times().size() - 1, s.get_n_rejections(), after.n_jac,
          after.n_rhs};
}

double relative_error(const State& a, const State& b) {
  double worst = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    worst = std::max(worst, std::abs(a[i] - b[i]) / std::abs(b[i]));
  }
  return worst;
}

// Kept across steps, the Jacobian is formed several times less often than
// RODAS forms it, at an error of the tolerance's order. Most of the saving is
// at a loose tolerance, where steps are long and few are rejected; at a tight
// one the third-order method takes more steps than RODAS's fourth, and more
// of them are rejected, each forming a new J.
void test_fewer_jacobians() {
  std::printf("Jacobians formed, against RODAS\n");
  const std::size_t age = ode::OdeControl().get_jacobian_max_age();
  for (const problem& p :
       {problem{"Robertson", robertson, robertson_jac, {1.0, 0.0, 0.0}, 40.0,
                1e-6},
        problem{"Van der Pol", vanderpol, vanderpol_jac, {2.0, 0.0}, 2.0,
                1.0}}) {
    const State reference = solve(p, Method::rodas, 1e-12, 1).y;
    for (const auto& [tol, fewer] :
         {std::pair{1e-2, 3.0}, std::pair{1e-4, 2.0}, std::pair{1e-6, 2.0}}) {
      const outcome r = solve(p, Method::rodas, tol, 1);
      const outcome w = solve(p, Method::ros34pw2, tol, age);
      const double err = relative_error(w.y, reference);
      std::printf("       (%s at %.0e: RODAS %zu J, %zu steps; ROS34PW2 %zu J, "
                  "%zu steps, %zu rejected, error %.1e)\n",
                  p.name, tol, r.jac, r.accepted, w.jac, w.accepted,
                  w.rejected, err);
      char label[96];
      std::snprintf(label, sizeof label,
                    "%s at %.0e: ROS34PW2 forms J at least %.0f times less "
                    "often than RODAS", p.name, tol, fewer);
      check(w.jac * fewer <= r.jac, label);
      check(err < 10.0 * tol, "  and its error is of the tolerance's order");
    }
  }
}

// When J is formed: at the first step, at a retry from a point J was not
// formed at, and when it has served `age` steps. So with age 1 it is formed
// once per accepted step and never for a retry from that step's start (J was
// formed there), and with age k it is formed at least once per k accepted
// steps and at most that plus once per rejection.
void test_when_formed() {
  std::printf("when the Jacobian is formed\n");
  const problem p{"Van der Pol", vanderpol, vanderpol_jac, {2.0, 0.0}, 2.0, 1.0};
  const outcome every = solve(p, Method::ros34pw2, 1e-4, 1);
  check(every.rejected > 0, "the run rejects steps (the test is not vacuous)");
  check(every.jac == every.accepted,
        "age 1: one J per accepted step, none for a retry");
  for (std::size_t age : {5, 20}) {
    const outcome w = solve(p, Method::ros34pw2, 1e-4, age);
    const std::size_t at_least = (w.accepted + age - 1) / age;
    check(w.jac >= at_least && w.jac <= at_least + w.rejected,
          "age " + std::to_string(age) +
            ": between one J per `age` steps and that plus one per rejection");
  }
  // 4 evaluations per attempt (three stages and the new point), and one at
  // the start; df/dt is not asked of an autonomous system.
  check(every.rhs == 1 + 4 * (every.accepted + every.rejected),
        "four right-hand sides per attempt");
}

// A state loaded with set_state starts a run as a new solver does: it keeps
// no Jacobian from the run before.
void test_set_state_starts_afresh() {
  std::printf("set_state\n");
  const problem p{"Van der Pol", vanderpol, vanderpol_jac, {2.0, 0.0}, 2.0, 1.0};
  CallbackSystem sys(p.rhs, p.y0, 0.0, p.jac, CallbackSystem::valid_type(),
                     true);
  ode::OdeControl control(1e-4, 1e-4, 1.0, 0.0, 1e-14, 2.0, 1e-6);
  ode::Solver<CallbackSystem> fresh(sys, control, Method::ros34pw2);
  ode::Solver<CallbackSystem> reused(sys, control, Method::ros34pw2);
  fresh.set_collect(false);
  reused.set_collect(false);
  reused.advance_adaptive({0.0, 1.0});
  reused.set_state(p.y0, 0.0);
  fresh.advance_adaptive({0.0, 2.0});
  reused.advance_adaptive({0.0, 2.0});
  check(fresh.state() == reused.state() && fresh.times() == reused.times(),
        "a run after set_state is a new solver's, to the bit");
}

// No sweep, no replay of a row, and no forward walk: the steps at which a run
// formed J are not recorded, so no walk could take the run's steps again.
void test_refusals() {
  std::printf("refusals\n");
  using Lorenz = LorenzSystem<double>;
  ode::Solver<Lorenz> s(Lorenz(10.0, 28.0, 8.0 / 3.0), ode::OdeControl(),
                        Method::ros34pw2);
  s.set_collect(false);
  s.set_keep_states(true);
  s.advance_adaptive({0.0, 0.5});
  check(refused([&] {
          ode::adjoint_batch lambda = ode::adjoint_batch::one({1.0, 0.0, 0.0});
          ode::adjoint_batch dp(1, 3);
          s.sweep_adjoint(lambda, dp);
        }),
        "a sweep is refused");
  check(refused([&] {
          const std::size_t last = s.recording().size() - 1;
          s.sweep_tangent<1>(ode::tangent_directions{{0}, {}}, {last});
        }),
        "and so is a forward walk");
}

}  // namespace

int main() {
  std::printf("odelia Rosenbrock-W stepper, standalone\n");
  test_order();
  test_fewer_jacobians();
  test_when_formed();
  test_set_state_starts_afresh();
  test_refusals();
  if (failures > 0) {
    std::printf("%d FAILED\n", failures);
    return 1;
  }
  std::printf("all passed\n");
  return 0;
}
