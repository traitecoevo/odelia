// Regression guard for the R-free solver core (#43).
//
// This translation unit includes every header in inst/include/odelia/ that is
// part of the reusable C++ core, and is compiled and linked with NO R and NO
// Rcpp anywhere -- see the Makefile alongside it. If someone adds an
// `#include <Rcpp.h>` or `<RcppCommon.h>` to one of these headers, or comes to
// rely on one arriving transitively, this stops building.
//
// It also catches the subtler version of the same bug: a header using a
// standard-library facility (`assert`, `std::string`, ...) that it never
// includes and only receives by accident from R's headers. That is exactly how
// spline.hpp came to use `assert` without <cassert>.
//
// The R interface headers -- solver_interface.hpp, rcpp_interface_helpers.hpp
// -- are deliberately NOT listed here. They are meant to depend on Rcpp.
//
// Note that this links src/Tape.cpp: the XAD tape runtime lives in exactly one
// object file, by design (see ARCHITECTURE.md), so a standalone consumer that
// instantiates Solver has to compile it too. That is a linking requirement, not
// an R one.

#include <odelia/ode_util.hpp>
#include <odelia/spline.hpp>
#include <odelia/interpolator.hpp>
#include <odelia/drivers.hpp>
#include <odelia/ode_control.hpp>
#include <odelia/ode_interface.hpp>
#include <odelia/ode_step.hpp>
#include <odelia/ode_step_rodas.hpp>
#include <odelia/ode_solver_internal.hpp>
#include <odelia/ode_solver.hpp>
#include <odelia/ode_fit.hpp>
#include <examples/lorenz_system.hpp>

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int failures = 0;

void check(bool ok, const std::string &what) {
  std::printf("  %s %s\n", ok ? "ok  " : "FAIL", what.c_str());
  if (!ok) {
    ++failures;
  }
}

// util::stop must throw rather than reach for R, and must carry its message.
void test_stop_throws() {
  bool threw = false;
  std::string msg;
  try {
    odelia::util::check_length(2, 3);
  } catch (const std::runtime_error &e) {
    threw = true;
    msg = e.what();
  }
  check(threw, "check_length throws std::runtime_error");
  check(msg.find("expected 3, received 2") != std::string::npos,
        "the message survives the throw");
}

// The interpolator is the part of the core most downstream consumers touch.
void test_interpolator() {
  odelia::interpolator::Interpolator in;
  in.init({0.0, 1.0, 2.0, 3.0}, {0.0, 1.0, 4.0, 9.0});
  check(std::abs(in.eval(2.0) - 4.0) < 1e-12, "interpolator hits its knots");

  bool threw = false;
  try {
    odelia::interpolator::Interpolator too_short;
    too_short.init({0.0, 1.0}, {0.0, 1.0});
  } catch (const std::runtime_error &) {
    threw = true;
  }
  check(threw, "interpolator rejects fewer than three points");
}

// Integrating a real system with no R session anywhere is the whole point.
// Tolerance-free check: tightening the controller must not move the answer.
void test_solver_runs() {
  using System = LorenzSystem<double>;
  const std::vector<double> times{0.0, 0.5, 1.0};

  System sys(10.0, 28.0, 8.0 / 3.0);
  odelia::ode::Solver<System> solver(sys, odelia::ode::OdeControl());
  solver.advance_adaptive(times);
  const auto loose = solver.state();

  odelia::ode::OdeControl tight;
  tight.set_tol_abs(1e-12);
  tight.set_tol_rel(1e-12);
  System sys_tight(10.0, 28.0, 8.0 / 3.0);
  odelia::ode::Solver<System> solver_tight(sys_tight, tight);
  solver_tight.advance_adaptive(times);
  const auto converged = solver_tight.state();

  check(solver.time() == 1.0, "solver arrives at the requested time");
  check(loose.size() == 3 && converged.size() == 3, "state has three elements");

  double worst = 0.0;
  for (size_t i = 0; i < loose.size(); ++i) {
    worst = std::max(worst, std::abs(loose[i] - converged[i]));
  }
  check(worst < 1e-4, "the two tolerances agree on the trajectory");
}

// A non-finite error estimate means the step left the model's valid domain, and
// must be rejected. It used to be *accepted*: NaN compares false against
// everything, so `rmax > 1.1` and `rmax < 0.5` were both false and control fell
// through to the branch that reports a successful step (odelia#52).
void test_control_rejects_nonfinite_error() {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  const std::vector<double> y{0.2, 0.2, 0.2};
  const std::vector<double> dydt{0.0, 0.0, 0.0};
  const double h = 1e-3;

  // tol_abs, tol_rel, a_y, a_dydt, h_min, h_max, h_init
  auto control = [] {
    return odelia::ode::OdeControl(1e-4, 1e-4, 1.0, 0.0, 1e-6, 5.0, 1e-6);
  };

  {
    // Position used to matter, which is why this went unnoticed. The reduction
    // was `rmax = std::max(r, rmax)`, and std::max(a, b) is `(a < b) ? b : a`:
    // with a = NaN it returns NaN, but on the *next* element a finite r returns
    // r, wiping the NaN. So a NaN only reached the branch if it survived to the
    // end of the loop -- last element, or all of them. In plant that is exactly
    // the case that bites: Patch chains the environment state *after* the
    // species, so the soil block sits in the trailing indices.
    odelia::ode::OdeControl c = control();
    const std::vector<double> yerr{1e-3, 1e-3, nan};   // NaN last: the real bug
    const double next = c.adjust_step_size(y.size(), 5, h, y, yerr, dydt);
    check(c.step_size_shrank(), "a trailing NaN error estimate rejects the step");
    check(next < h, "a trailing NaN error estimate shrinks the step");
  }
  {
    // A NaN anywhere must reject, not just where the reduction happened to
    // preserve it. Passed before the fix too, but for the wrong reason -- the
    // trailing finite element rejected on its own magnitude.
    odelia::ode::OdeControl c = control();
    const std::vector<double> yerr{1e-3, nan, 1e-3};
    c.adjust_step_size(y.size(), 5, h, y, yerr, dydt);
    check(c.step_size_shrank(), "an interior NaN error estimate rejects the step");
  }
  {
    // The sharpest form: a NaN alongside errors that would otherwise be
    // comfortably *accepted*. Before the fix the NaN was wiped and the step
    // grew.
    odelia::ode::OdeControl c = control();
    const std::vector<double> yerr{nan, 1e-12, 1e-12};
    const double next = c.adjust_step_size(y.size(), 5, h, y, yerr, dydt);
    check(c.step_size_shrank() && next < h,
          "a NaN masked by small finite errors still rejects");
  }
  {
    odelia::ode::OdeControl c = control();
    const std::vector<double> yerr{1e-3, inf, 1e-3};
    c.adjust_step_size(y.size(), 5, h, y, yerr, dydt);
    check(c.step_size_shrank(), "an Inf error estimate rejects the step");
  }
  {
    // A NaN in the *state* poisons errlevel, so it arrives as a NaN ratio too.
    odelia::ode::OdeControl c = control();
    const std::vector<double> y_bad{0.2, nan, 0.2};
    const std::vector<double> yerr{1e-3, 1e-3, 1e-3};
    c.adjust_step_size(y_bad.size(), 5, h, y_bad, yerr, dydt);
    check(c.step_size_shrank(), "a NaN state rejects the step");
  }
  {
    // Already at the floor: cannot decrease, but must still report the shrink
    // so the caller raises rather than committing the non-finite state.
    odelia::ode::OdeControl c = control();
    const std::vector<double> yerr{nan, nan, nan};
    const double next = c.adjust_step_size(y.size(), 5, 1e-6, y, yerr, dydt);
    check(c.step_size_shrank() && next == 1e-6,
          "at step_size_min a NaN still reports a shrink");
  }

  // Regressions: finite behaviour must be untouched.
  {
    odelia::ode::OdeControl c = control();
    const std::vector<double> yerr{1e-3, 1e-3, 1e-3};
    const double next = c.adjust_step_size(y.size(), 5, h, y, yerr, dydt);
    check(c.step_size_shrank() && next < h, "a large finite error still rejects");
  }
  {
    odelia::ode::OdeControl c = control();
    const std::vector<double> yerr{1e-12, 1e-12, 1e-12};
    const double next = c.adjust_step_size(y.size(), 5, h, y, yerr, dydt);
    check(!c.step_size_shrank() && next > h, "a small finite error still grows");
  }
}

// End to end: a system whose derivatives go non-finite outside a bounded range
// must fail loudly rather than integrate on with a poisoned state.
namespace {
struct BoundedSystem {
  using value_type = double;
  double y = 0.5, dydt = 0.0, time = 0.0;

  size_t ode_size() const { return 1; }
  double ode_time() const { return time; }

  template <typename Iterator> Iterator set_ode_state(Iterator it, double t) {
    y = *it++;
    time = t;
    // Valid only on [0, 1]; outside it the model has nothing to say.
    dydt = (y >= 0.0 && y <= 1.0)
               ? 50.0
               : std::numeric_limits<double>::quiet_NaN();
    return it;
  }
  template <typename Iterator> Iterator ode_state(Iterator it) const {
    *it++ = y;
    return it;
  }
  template <typename Iterator> Iterator ode_rates(Iterator it) const {
    *it++ = dydt;
    return it;
  }
  template <typename Iterator> Iterator ode_aux(Iterator it) const { return it; }
};
} // namespace

void test_solver_refuses_nonfinite_state() {
  BoundedSystem sys;
  odelia::ode::OdeControl control(1e-4, 1e-4, 1.0, 0.0, 1e-6, 5.0, 1.0);
  odelia::ode::Solver<BoundedSystem> solver(sys, control);

  bool threw = false;
  std::string msg;
  try {
    // One step of h = 1 at dydt = 50 lands far outside [0, 1].
    solver.advance_adaptive(std::vector<double>{0.0, 1.0});
  } catch (const std::runtime_error &e) {
    threw = true;
    msg = e.what();
  }

  const auto state = solver.state();
  const bool finite = state.empty() || std::isfinite(state[0]);
  check(threw || finite, "the solver never commits a non-finite state");
  if (threw) {
    std::printf("       (raised: %s)\n", msg.c_str());
  }
}

// --- Invariant-aware step rejection (#55) ----------------------------------
//
// The logistic flow dy/dt = k*y*(1-y) keeps the exact solution inside (0, 1) for
// any y0 in (0, 1), but a finite explicit step from y ~ 1/2 with k large lands
// well outside it. So the *discretisation* violates a bound the model itself
// respects, which is precisely the case this feature exists for -- and note that
// the rate stays perfectly finite out there, so nothing about it is visible to
// the non-finite check added for #52.
//
// OnLeave picks how the system reports having left [0, 1].
namespace {

enum class OnLeave { nothing, throw_domain, throw_bug };

template <OnLeave on_leave>
struct Logistic {
  using value_type = double;
  static constexpr double k = 50.0;
  double y = 0.5, dydt = 0.0, time = 0.0;

  size_t ode_size() const { return 1; }
  double ode_time() const { return time; }

  template <typename Iterator> Iterator set_ode_state(Iterator it, double t) {
    y = *it++;
    time = t;
    if (y < 0.0 || y > 1.0) {
      if (on_leave == OnLeave::throw_domain) {
        odelia::util::stop_domain("y = " + odelia::util::to_string(y) +
                                  " is outside [0, 1]");
      } else if (on_leave == OnLeave::throw_bug) {
        odelia::util::stop("bug-shaped failure, must not be absorbed");
      }
    }
    dydt = k * y * (1.0 - y);
    return it;
  }
  template <typename Iterator> Iterator ode_state(Iterator it) const {
    *it++ = y;
    return it;
  }
  template <typename Iterator> Iterator ode_rates(Iterator it) const {
    *it++ = dydt;
    return it;
  }
  template <typename Iterator> Iterator ode_aux(Iterator it) const { return it; }
};

// Same dynamics, but declaring the domain. Inherited rather than switched on a
// template parameter so that the silent twin genuinely lacks the method and
// has_state_check<> resolves to false for it -- an `if constexpr` inside one
// struct would still leave the member there for the trait to find.
struct LogisticChecked : Logistic<OnLeave::nothing> {
  static int refusals;
  bool ode_state_valid(const std::vector<double>& state) const {
    const bool ok = state[0] >= 0.0 && state[0] <= 1.0;
    if (!ok) {
      ++refusals;
    }
    return ok;
  }
};
int LogisticChecked::refusals = 0;

// Loose enough that the error estimate alone is content to accept a step that
// leaves the domain, with a first step large enough to do so.
odelia::ode::OdeControl loose_control() {
  return odelia::ode::OdeControl(1e-2, 1e-2, 1.0, 0.0, 1e-8, 10.0, 1.0);
}

} // namespace

// A declared domain must be enforced on the committed state.
void test_predicate_rejects_out_of_domain_step() {
  check(odelia::ode::has_state_check<LogisticChecked>::value,
        "has_state_check finds a declared ode_state_valid");
  check(!odelia::ode::has_state_check<Logistic<OnLeave::nothing>>::value,
        "and does not invent one that is absent");

  Logistic<OnLeave::nothing> unguarded;
  odelia::ode::Solver<Logistic<OnLeave::nothing>> s0(unguarded, loose_control());
  s0.advance_adaptive(std::vector<double>{0.0, 1.0});

  LogisticChecked::refusals = 0;
  LogisticChecked guarded;
  odelia::ode::Solver<LogisticChecked> s1(guarded, loose_control());
  s1.advance_adaptive(std::vector<double>{0.0, 1.0});
  const double y = s1.state()[0];

  check(LogisticChecked::refusals > 0,
        "the predicate actually refused at least one step (test is not vacuous)");
  check(y >= 0.0 && y <= 1.0, "the committed state stays inside [0, 1]");
  check(s1.time() == 1.0, "and the solve still reaches the requested time");
  std::printf("       (%d refusal(s); unguarded twin finished at y = %g)\n",
              LogisticChecked::refusals, s0.state()[0]);
}

// The usual way a model reports an impossible state is to throw. That must cost
// the step, not the solve.
void test_domain_error_becomes_a_rejection() {
  Logistic<OnLeave::throw_domain> sys;
  odelia::ode::Solver<Logistic<OnLeave::throw_domain>> s(sys, loose_control());

  bool threw = false;
  std::string msg;
  try {
    s.advance_adaptive(std::vector<double>{0.0, 1.0});
  } catch (const std::runtime_error& e) {
    threw = true;
    msg = e.what();
  }

  check(!threw, "a DomainError from a stage is a rejection, not a fatal");
  if (threw) {
    std::printf("       (raised: %s)\n", msg.c_str());
  } else {
    const double y = s.state()[0];
    check(y >= 0.0 && y <= 1.0, "and the solve finishes inside the domain");
    check(s.time() == 1.0, "and reaches the requested time");
  }
}

// The other half of that bargain: only DomainError is absorbed. A plain
// util::stop() is how the core reports a bug, and turning one into step-shrinking
// would hide it behind an accuracy complaint.
void test_non_domain_throw_is_not_absorbed() {
  Logistic<OnLeave::throw_bug> sys;
  odelia::ode::Solver<Logistic<OnLeave::throw_bug>> s(sys, loose_control());

  bool threw = false;
  std::string msg;
  try {
    s.advance_adaptive(std::vector<double>{0.0, 1.0});
  } catch (const std::runtime_error& e) {
    threw = true;
    msg = e.what();
  }

  check(threw, "a util::stop() from a stage still ends the solve");
  check(msg.find("bug-shaped") != std::string::npos,
        "and arrives with its own message, not an accuracy complaint");
}

namespace {
// dy/dt = 1 under a ceiling of y <= 1, started at 0.9 and asked to advance a full
// unit of time: the *exact* flow leaves the domain, so no step size helps. Linear,
// so the RKCK error estimate is essentially zero and the accuracy controller never
// interferes -- the predicate is the only thing that can reject.
struct RampToCeiling {
  using value_type = double;
  double y = 0.9, dydt = 1.0, time = 0.0;

  size_t ode_size() const { return 1; }
  double ode_time() const { return time; }

  template <typename Iterator> Iterator set_ode_state(Iterator it, double t) {
    y = *it++;
    time = t;
    dydt = 1.0;
    return it;
  }
  template <typename Iterator> Iterator ode_state(Iterator it) const {
    *it++ = y;
    return it;
  }
  template <typename Iterator> Iterator ode_rates(Iterator it) const {
    *it++ = dydt;
    return it;
  }
  template <typename Iterator> Iterator ode_aux(Iterator it) const { return it; }

  bool ode_state_valid(const std::vector<double>& state) const {
    return state[0] <= 1.0;
  }
};
} // namespace

// Rejection cannot rescue a model whose exact flow leaves its domain -- it can
// only detect it. It must then fail saying so, rather than blaming accuracy.
void test_unreachable_domain_fails_with_a_reason() {
  RampToCeiling sys;
  odelia::ode::OdeControl control(1e-6, 1e-6, 1.0, 0.0, 1e-8, 1.0, 0.5);
  odelia::ode::Solver<RampToCeiling> s(sys, control);

  bool threw = false;
  std::string msg;
  try {
    s.advance_adaptive(std::vector<double>{0.0, 1.0});
  } catch (const std::runtime_error& e) {
    threw = true;
    msg = e.what();
  }

  check(threw, "a domain the exact flow leaves cannot be integrated");
  check(msg.find("invalid state") != std::string::npos &&
            msg.find("ode_state_valid") != std::string::npos,
        "and the failure names the reason rather than blaming accuracy");
  if (threw) {
    std::printf("       (raised: %s)\n", msg.c_str());
  }
}

} // namespace

int main() {
  std::printf("odelia solver core, standalone (no R, no Rcpp)\n");
  test_stop_throws();
  test_interpolator();
  test_solver_runs();
  test_control_rejects_nonfinite_error();
  test_solver_refuses_nonfinite_state();
  test_predicate_rejects_out_of_domain_step();
  test_domain_error_becomes_a_rejection();
  test_non_domain_throw_is_not_absorbed();
  test_unreachable_domain_fails_with_a_reason();
  if (failures == 0) {
    std::printf("all checks passed\n");
    return 0;
  }
  std::printf("%d failure(s)\n", failures);
  return 1;
}
