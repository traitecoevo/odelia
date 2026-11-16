// Event location (#92): Solver::advance_until, and the two pieces it is built
// from, SolverInternal::locate and retract_step. What it finds, that the state
// it stops on is the one a run asked to land on the event would hold, and that
// the search on the interpolant allocates nothing.

#include <odelia/ode_solver.hpp>

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <numbers>
#include <optional>
#include <string>
#include <vector>

// Every allocation in this program is counted, so a test can require that a
// search makes none.
static std::size_t allocations = 0;

void* operator new(std::size_t n) {
  ++allocations;
  if (void* p = std::malloc(n == 0 ? 1 : n)) {
    return p;
  }
  throw std::bad_alloc();
}
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }

namespace {

namespace ode = odelia::ode;
namespace root = odelia::root;
using ode::Method;

int failures = 0;

void check(bool ok, const std::string& what) {
  std::printf("  %s %s\n", ok ? "ok  " : "FAIL", what.c_str());
  if (!ok) {
    ++failures;
  }
}

// Equal element for element, a NaN equal to a NaN: the first step size, which
// no step took, is NaN in every run.
bool same(const std::vector<double>& a, const std::vector<double>& b) {
  if (a.size() != b.size()) {
    return false;
  }
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (!(a[i] == b[i] || (std::isnan(a[i]) && std::isnan(b[i])))) {
      return false;
    }
  }
  return true;
}

const char* name(Method m) {
  return m == Method::rkck ? "rkck" : m == Method::dopri ? "dopri" : "rodas";
}

// x' = v, v' = -x from (1, 0): x = cos t, which crosses zero at pi/2 + k pi.
struct Oscillator {
  using value_type = double;
  double x = 1.0, v = 0.0, time = 0.0;

  std::size_t ode_size() const { return 2; }
  double ode_time() const { return time; }
  template <typename It> It set_ode_state(It it, double t) {
    x = *it++;
    v = *it++;
    time = t;
    return it;
  }
  template <typename It> It ode_state(It it) const {
    *it++ = x;
    *it++ = v;
    return it;
  }
  template <typename It> It ode_rates(It it) const {
    *it++ = v;
    *it++ = -x;
    return it;
  }
  template <typename It> It ode_aux(It it) const { return it; }
  // For the Rosenbrock stepper.
  void ode_jacobian(const std::vector<double>&, double,
                    const std::vector<double>&, std::vector<double>& J) {
    J.assign({0.0, 1.0, -1.0, 0.0});
  }
};

// y' = y from 1: y reaches L at log(L).
struct Growth {
  using value_type = double;
  double y = 1.0, time = 0.0;

  std::size_t ode_size() const { return 1; }
  double ode_time() const { return time; }
  template <typename It> It set_ode_state(It it, double t) {
    y = *it++;
    time = t;
    return it;
  }
  template <typename It> It ode_state(It it) const {
    *it++ = y;
    return it;
  }
  template <typename It> It ode_rates(It it) const {
    *it++ = y;
    return it;
  }
  template <typename It> It ode_aux(It it) const { return it; }
  void ode_jacobian(const std::vector<double>&, double,
                    const std::vector<double>&, std::vector<double>& J) {
    J.assign({1.0});
  }
};

ode::OdeControl control(ode::Controller rule) {
  ode::OdeControl c(1e-8, 1e-8, 1.0, 0.0, 1e-8, 10.0, 1e-6);
  c.set_controller(rule);
  return c;
}

const root::tolerance tight{1e-12, 1e-12};

const auto x_is_zero = [](double, const std::vector<double>& y) { return y[0]; };

void test_finds_the_crossing(Method m) {
  ode::Solver<Oscillator> s(Oscillator{}, control(ode::Controller::hairer), m);
  s.set_collect(false);
  const std::optional<double> t = s.advance_until(x_is_zero, 10.0, tight);
  const double tol = m == Method::rodas ? 1e-6 : 1e-7;
  check(t.has_value() && std::abs(*t - std::numbers::pi / 2) < tol,
        std::string(name(m)) + ": x = cos t first reaches zero at pi/2");
  check(t.has_value() && s.time() == *t,
        "and the solver stops at the located time exactly");
  check(std::abs(s.state()[0]) < tol && std::abs(s.state()[1] + 1.0) < tol,
        "holding the state there, x = 0 and v = -1");
}

// The state at the event is the one a run asked to land on that time holds:
// the same steps to the start of the step that crossed, then steps to the
// event under the same control. Run B takes A's steps by copying the solver
// before each, and from the copy before the step that passed the event, steps
// to it. Everything after is compared too.
void test_is_a_landing_time(Method m, ode::Controller rule) {
  ode::Solver<Oscillator> a(Oscillator{}, control(rule), m);
  a.set_collect(false);
  const std::optional<double> t = a.advance_until(x_is_zero, 10.0, tight);
  if (!t) {
    check(false, "no event found");
    return;
  }

  ode::Solver<Oscillator> b(Oscillator{}, control(rule), m);
  b.set_collect(false);
  while (true) {
    ode::Solver<Oscillator> before = b;
    b.step(10.0);
    if (b.time() >= *t) {
      b = before;
      break;
    }
  }
  while (b.time() < *t) {
    b.step(*t);
  }
  const std::string label = std::string(name(m)) + ", " +
      (rule == ode::Controller::hairer ? "hairer" : "gsl") + ": ";
  check(a.time() == b.time() && a.state() == b.state() &&
            a.times() == b.times() && same(a.step_sizes(), b.step_sizes()) &&
            a.get_step_size() == b.get_step_size(),
        label + "the event's state, times and next step are a landing run's");
  a.advance_adaptive({a.time(), 10.0});
  b.advance_adaptive({b.time(), 10.0});
  check(a.state() == b.state() && a.times() == b.times() &&
            same(a.step_sizes(), b.step_sizes()),
        label + "and so is every step after it");
}

void test_no_event_is_an_ordinary_run(Method m) {
  ode::Solver<Growth> a(Growth{}, control(ode::Controller::hairer), m);
  ode::Solver<Growth> b(Growth{}, control(ode::Controller::hairer), m);
  const auto past_100 = [](double, const std::vector<double>& y) {
    return y[0] - 100.0;
  };
  const std::optional<double> t = a.advance_until(past_100, 2.0, tight);
  b.advance_adaptive({0.0, 2.0});
  check(!t && a.time() == 2.0,
        std::string(name(m)) + ": with no event before time_max it stops there");
  check(a.state() == b.state() && a.times() == b.times(),
        "exactly as advance_adaptive() does");
}

void test_edges() {
  {
    ode::Solver<Growth> s(Growth{}, control(ode::Controller::hairer));
    const auto at_one = [](double, const std::vector<double>& y) {
      return y[0] - 1.0;
    };
    const std::optional<double> t = s.advance_until(at_one, 2.0, tight);
    check(t && *t == 0.0 && s.times().size() == 1,
          "g = 0 at the start is an event there, and no step is taken");
  }
  {
    ode::Solver<Growth> s(Growth{}, control(ode::Controller::hairer));
    const auto at_end = [](double t, const std::vector<double>&) {
      return t - 2.0;
    };
    const std::optional<double> t = s.advance_until(at_end, 2.0, tight);
    check(t && *t == 2.0 && s.time() == 2.0,
          "g reaching zero exactly at a step's end is an event at that end");
  }
  {
    ode::Solver<Growth> s(Growth{}, control(ode::Controller::hairer));
    bool threw = false;
    try {
      s.advance_until(x_is_zero, -1.0, tight);
    } catch (const std::runtime_error&) {
      threw = true;
    }
    check(threw, "a time_max before the current time is refused");
  }
}

// A sequence of levels, each from where the last stopped: what a run to a
// series of sizes asks.
void test_a_sequence_of_levels(Method m) {
  ode::Solver<Growth> s(Growth{}, control(ode::Controller::hairer), m);
  s.set_collect(false);
  double worst = 0.0;
  bool all = true;
  for (double level = 2.0; level <= 10.0; level += 1.0) {
    const auto reaches = [level](double, const std::vector<double>& y) {
      return y[0] - level;
    };
    const std::optional<double> t = s.advance_until(reaches, 100.0, tight);
    all = all && t.has_value();
    if (t) {
      worst = std::max(worst, std::abs(*t - std::log(level)));
    }
  }
  char worst_text[32];
  std::snprintf(worst_text, sizeof worst_text, "%.1e", worst);
  check(all && worst < 1e-6,
        std::string(name(m)) + ": y = e^t reaches 2, ..., 10 at their logs, "
        "worst error " + worst_text);
}

void test_locate_allocates_nothing() {
  Oscillator sys;
  ode::SolverInternal<Oscillator> s(sys, control(ode::Controller::hairer),
                                    Method::dopri);
  double g_prev = 1.0;
  while (true) {
    s.step(sys, 10.0);
    if (s.state_ref()[0] < 0.0) {
      break;
    }
    g_prev = s.state_ref()[0];
  }
  const double g_end = s.state_ref()[0];
  s.locate(x_is_zero, g_prev, g_end, tight);  // sizes the scratch
  const std::size_t before = allocations;
  double sum = 0.0;
  for (int i = 0; i < 100; ++i) {
    sum += s.locate(x_is_zero, g_prev, g_end, tight);
  }
  const std::size_t made = allocations - before;
  check(made == 0 && std::isfinite(sum),
        "100 locations on the interpolant made " + std::to_string(made) +
            " allocations");
}

void test_retract_refusals() {
  Growth sys;
  ode::SolverInternal<Growth> s(sys, control(ode::Controller::hairer),
                                Method::dopri);
  const auto refused = [&](auto&& call) {
    try {
      call();
    } catch (const std::runtime_error&) {
      return true;
    }
    return false;
  };
  check(refused([&] { s.retract_step(sys); }),
        "retract_step() with nothing held is refused");
  s.hold_step_start();
  s.step(sys, 1.0);
  s.step(sys, 1.0);
  check(refused([&] { s.retract_step(sys); }),
        "and so is retracting two steps");
  s.hold_step_start();
  s.step_to(sys, 1.5);
  check(refused([&] { s.retract_step(sys); }),
        "and a pinned step, which leaves no interpolant");
}

}  // namespace

int main() {
  std::printf("odelia event location, standalone\n");
  for (Method m : {Method::rkck, Method::dopri, Method::rodas}) {
    test_finds_the_crossing(m);
  }
  for (Method m : {Method::rkck, Method::dopri, Method::rodas}) {
    for (ode::Controller rule : {ode::Controller::hairer, ode::Controller::gsl}) {
      test_is_a_landing_time(m, rule);
    }
  }
  for (Method m : {Method::rkck, Method::dopri, Method::rodas}) {
    test_no_event_is_an_ordinary_run(m);
    test_a_sequence_of_levels(m);
  }
  test_edges();
  test_locate_allocates_nothing();
  test_retract_refusals();
  if (failures > 0) {
    std::printf("%d FAILED\n", failures);
    return 1;
  }
  std::printf("all passed\n");
  return 0;
}
