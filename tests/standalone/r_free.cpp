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
#include <odelia/ode_step_dopri.hpp>
#include <odelia/ode_solver_internal.hpp>
#include <odelia/ode_solver.hpp>
#include <odelia/ode_fit.hpp>
#include <odelia/ode_callback_system.hpp>
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

// --- The same bargain on the pinned path (plant#642) ------------------------
//
// advance_fixed() steps exactly to a caller-supplied set of times, which is how
// a replay reproduces a trajectory recorded earlier. Its endpoints cannot move,
// so #55's answer to an invalid step -- take a smaller one -- has to become
// "take several smaller ones to the same endpoint". Until it did, advance_fixed
// called the stepper bare and the first refusal killed the solve.

namespace {
// The grid below is coarse enough that a single RKCK step from y = 0.5 leaves
// [0, 1], so every one of these tests exercises the subdivision rather than
// merely passing through it.
std::vector<double> coarse_grid() {
  return std::vector<double>{0.0, 0.25, 0.5, 0.75, 1.0};
}
} // namespace

void test_pinned_step_domain_error_is_a_rejection() {
  Logistic<OnLeave::throw_domain> sys;
  odelia::ode::Solver<Logistic<OnLeave::throw_domain>> s(sys, loose_control());

  bool threw = false;
  std::string msg;
  try {
    s.advance_fixed(coarse_grid());
  } catch (const std::runtime_error& e) {
    threw = true;
    msg = e.what();
  }

  check(!threw, "a DomainError under advance_fixed costs the step, not the solve");
  if (threw) {
    std::printf("       (raised: %s)\n", msg.c_str());
    return;
  }
  const double y = s.state()[0];
  check(y >= 0.0 && y <= 1.0, "and the solve finishes inside the domain");
  // The whole point of the pinned path: subdividing must not move the endpoint,
  // because the caller matches these times against its own record.
  check(s.time() == coarse_grid().back(),
        "and lands exactly on the last requested time");
}

// The predicate is the other way #55 lets a system refuse a state, and it must
// be enforced here too -- otherwise it would go silently unchecked whenever the
// integration happened to be pinned.
void test_pinned_step_predicate_is_enforced() {
  LogisticChecked sys;
  odelia::ode::Solver<LogisticChecked> s(sys, loose_control());

  const int refusals_before = LogisticChecked::refusals;
  bool threw = false;
  try {
    s.advance_fixed(coarse_grid());
  } catch (const std::runtime_error&) {
    threw = true;
  }

  check(!threw, "a refused state under advance_fixed is a rejection, not a fatal");
  check(LogisticChecked::refusals > refusals_before,
        "and ode_state_valid() was actually consulted on the pinned path");
  if (!threw) {
    const double y = s.state()[0];
    check(y >= 0.0 && y <= 1.0, "the committed state stays inside [0, 1]");
    check(s.time() == coarse_grid().back(),
          "and the endpoint is still hit exactly");
  }
}

// Only DomainError is absorbed here, exactly as on the adaptive path.
void test_pinned_step_non_domain_throw_is_not_absorbed() {
  Logistic<OnLeave::throw_bug> sys;
  odelia::ode::Solver<Logistic<OnLeave::throw_bug>> s(sys, loose_control());

  bool threw = false;
  std::string msg;
  try {
    s.advance_fixed(coarse_grid());
  } catch (const std::runtime_error& e) {
    threw = true;
    msg = e.what();
  }

  check(threw, "a util::stop() from a stage still ends a pinned solve");
  check(msg.find("bug-shaped") != std::string::npos,
        "and arrives with its own message, not a subdivision complaint");
}

// Subdivision detects an unreachable domain; it cannot rescue one. Saying so is
// the difference between a diagnosis and an infinite loop.
void test_pinned_step_unreachable_domain_fails_with_a_reason() {
  RampToCeiling sys;
  odelia::ode::OdeControl control(1e-6, 1e-6, 1.0, 0.0, 1e-8, 1.0, 0.5);
  odelia::ode::Solver<RampToCeiling> s(sys, control);

  bool threw = false;
  std::string msg;
  try {
    s.advance_fixed(std::vector<double>{0.0, 0.5, 1.0});
  } catch (const std::runtime_error& e) {
    threw = true;
    msg = e.what();
  }

  check(threw, "a domain the exact flow leaves cannot be integrated pinned either");
  check(msg.find("invalid state") != std::string::npos &&
            msg.find("ode_state_valid") != std::string::npos,
        "and the failure names the reason and the sub-step it gave up at");
  if (threw) {
    std::printf("       (raised: %s)\n", msg.c_str());
  }
}

namespace {
// Counts derivative evaluations, to show that a system raising no objection is
// stepped exactly as it was before subdivision existed.
struct CountingLinear {
  using value_type = double;
  static int derivs;
  double y = 1.0, dydt = 1.0, time = 0.0;

  size_t ode_size() const { return 1; }
  double ode_time() const { return time; }

  template <typename Iterator> Iterator set_ode_state(Iterator it, double t) {
    y = *it++;
    time = t;
    dydt = 1.0;
    ++derivs;
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
int CountingLinear::derivs = 0;
} // namespace

// The cost of the new machinery on the ordinary path must be nothing: one RKCK
// step per interval, the same six stage evaluations, the same endpoints.
void test_pinned_step_unchanged_when_nothing_objects() {
  CountingLinear sys;
  odelia::ode::Solver<CountingLinear> s(sys, loose_control());

  const std::vector<double> times{0.0, 0.25, 0.5, 0.75, 1.0};
  CountingLinear::derivs = 0;
  s.advance_fixed(times);

  // Four intervals, six stage evaluations each. setup_dydt_in() adds none: the
  // solver starts with a clean dydt_in and FSAL keeps it clean thereafter.
  check(CountingLinear::derivs == 24,
        "an unobjecting system still takes exactly one step per interval");
  check(s.time() == 1.0, "and reaches the requested time");
  // dy/dt = 1 from y = 1 over one unit of time, integrated exactly by RKCK.
  check(std::abs(s.state()[0] - 2.0) < 1e-12, "with the exact trajectory");
  std::printf("       (%d derivative evaluation(s) over %zu intervals)\n",
              CountingLinear::derivs, times.size() - 1);
}


// --- A right-hand side handed in at run time (#62) ---------------------------
//
// CallbackSystem wraps a std::function as a System, so an R closure (or a
// Python callable, or a lambda as here) can be stepped by the same solver as a
// compiled system. The implicit stepper needs a Jacobian, which such a system
// supplies through the ode_jacobian() hook -- analytic or by finite differences
// -- rather than by AD on a rebound twin.

namespace {

using odelia::ode::CallbackSystem;
using State = std::vector<double>;

const double SIGMA = 10.0, RHO = 28.0, BETA = 8.0 / 3.0;

void lorenz_rhs(double, const State& y, State& dydt) {
  dydt[0] = SIGMA * (y[1] - y[0]);
  dydt[1] = y[0] * (RHO - y[2]) - y[1];
  dydt[2] = y[0] * y[1] - BETA * y[2];
}

void lorenz_jac(double, const State& y, const State&, State& J) {
  // row-major, J[row * 3 + col] = d f_row / d y_col
  J[0 * 3 + 0] = -SIGMA;    J[0 * 3 + 1] = SIGMA;  J[0 * 3 + 2] = 0.0;
  J[1 * 3 + 0] = RHO - y[2]; J[1 * 3 + 1] = -1.0;  J[1 * 3 + 2] = -y[0];
  J[2 * 3 + 0] = y[1];       J[2 * 3 + 1] = y[0];  J[2 * 3 + 2] = -BETA;
}

// A compiled system with its own Jacobian hook and no rebind(): the hook alone
// must open the implicit stepper to it.
struct HookedLorenz {
  using value_type = double;
  State y{1.0, 1.0, 1.0}, dydt{0.0, 0.0, 0.0};
  double time = 0.0;

  size_t ode_size() const { return 3; }
  double ode_time() const { return time; }
  template <typename Iterator> Iterator set_ode_state(Iterator it, double t) {
    for (size_t i = 0; i < 3; ++i) y[i] = *it++;
    time = t;
    lorenz_rhs(t, y, dydt);
    return it;
  }
  template <typename Iterator> Iterator ode_state(Iterator it) const {
    for (size_t i = 0; i < 3; ++i) *it++ = y[i];
    return it;
  }
  template <typename Iterator> Iterator ode_rates(Iterator it) const {
    for (size_t i = 0; i < 3; ++i) *it++ = dydt[i];
    return it;
  }
  template <typename Iterator> Iterator ode_aux(Iterator it) const { return it; }
  void ode_jacobian(const State& y_, double t, const State& f, State& J) const {
    lorenz_jac(t, y_, f, J);
  }
};

odelia::ode::OdeControl tight_control() {
  return odelia::ode::OdeControl(1e-10, 1e-10, 1.0, 0.0, 1e-12, 10.0, 1e-6);
}

double max_abs_diff(const State& a, const State& b) {
  double worst = 0.0;
  for (size_t i = 0; i < a.size(); ++i) {
    worst = std::max(worst, std::abs(a[i] - b[i]));
  }
  return worst;
}

} // namespace

void test_jacobian_hook_is_detected() {
  using odelia::ode::has_jacobian;
  using odelia::ode::has_rebind;
  using odelia::ode::Jacobian;
  using odelia::ode::RodasStep;

  check(has_jacobian<HookedLorenz>::value, "has_jacobian finds a declared ode_jacobian");
  check(!has_rebind<HookedLorenz>::value, "on a system with no rebind()");
  check(Jacobian<HookedLorenz>::supported && !Jacobian<HookedLorenz>::ad_supported,
        "so the Jacobian is supported through the hook, not AD");
  check(RodasStep<HookedLorenz>::supported, "and RODAS is open to it");

  check(has_jacobian<CallbackSystem>::value, "CallbackSystem declares the hook");
  check(RodasStep<CallbackSystem>::supported, "and RODAS is open to it too");

  check(!has_jacobian<LorenzSystem<double>>::value &&
            Jacobian<LorenzSystem<double>>::ad_supported,
        "the compiled Lorenz keeps its AD route (no hook, has rebind)");
  check(!Jacobian<Logistic<OnLeave::nothing>>::supported,
        "a system with neither hook nor rebind is not supported");

  check(odelia::ode::has_autonomous<CallbackSystem>::value,
        "has_autonomous finds CallbackSystem's declaration");
  check(!odelia::ode::has_autonomous<HookedLorenz>::value,
        "and does not invent one that is absent");
}

// The implicit stepper on a callback system must give the same trajectory as it
// does on the compiled system with its AD Jacobian, whether the hook is analytic
// or finite differences, and the same as the explicit stepper.
void test_callback_rodas_matches_compiled() {
  using odelia::ode::Method;
  const State y0{1.0, 1.0, 1.0};
  const std::vector<double> times{0.0, 1.0, 2.0};

  LorenzSystem<double> compiled(SIGMA, RHO, BETA);
  odelia::ode::Solver<LorenzSystem<double>> ref(compiled, tight_control(), Method::rodas);
  ref.advance_adaptive(times);
  const State reference = ref.state();

  CallbackSystem explicit_sys(lorenz_rhs, y0, 0.0);
  odelia::ode::Solver<CallbackSystem> rk(explicit_sys, tight_control(), Method::rkck);
  rk.advance_adaptive(times);
  check(max_abs_diff(rk.state(), reference) < 1e-5,
        "RKCK on a callback system agrees with compiled RODAS");

  CallbackSystem analytic(lorenz_rhs, y0, 0.0, lorenz_jac, CallbackSystem::valid_type(), true);
  odelia::ode::Solver<CallbackSystem> ra(analytic, tight_control(), Method::rodas);
  ra.advance_adaptive(times);
  check(max_abs_diff(ra.state(), reference) < 1e-5,
        "RODAS with an analytic Jacobian callback agrees");

  CallbackSystem fd(lorenz_rhs, y0, 0.0, CallbackSystem::jac_type(), CallbackSystem::valid_type(), true);
  odelia::ode::Solver<CallbackSystem> rf(fd, tight_control(), Method::rodas);
  rf.advance_adaptive(times);
  check(max_abs_diff(rf.state(), reference) < 1e-5,
        "RODAS with a finite-difference Jacobian agrees");

  HookedLorenz hooked;
  odelia::ode::Solver<HookedLorenz> rh(hooked, tight_control(), Method::rodas);
  rh.advance_adaptive(times);
  check(max_abs_diff(rh.state(), reference) < 1e-5,
        "RODAS through a compiled system's own hook agrees");

  // fd_jacobian itself, against the analytic matrix at a point.
  CallbackSystem probe(lorenz_rhs, y0, 0.0);
  const State y{1.5, -2.0, 20.0};
  State f(3), J_fd, J_an(9);
  odelia::ode::derivs(probe, y, f, 0.0);
  odelia::ode::fd_jacobian(probe, y, 0.0, f, J_fd);
  lorenz_jac(0.0, y, f, J_an);
  check(max_abs_diff(J_fd, J_an) < 1e-4, "fd_jacobian matches the analytic Jacobian");
}

// What a step costs in right-hand-side evaluations is the whole cost of a step
// for a consumer whose right-hand side is an equilibrium solve. The budget:
// one evaluation when the solver seeds itself, then per attempt five stages
// plus the derivative at the new point (6, whether the attempt is accepted or
// rejected -- the estimate that rejects it needs them), plus once per accepted
// step one for df/dt unless the system says it is autonomous, with the
// Jacobian (and df/dt) formed once per accepted step and kept across a retry.
void test_callback_call_budget() {
  using odelia::ode::Method;
  const State y0{1.0, 1.0, 1.0};
  const std::vector<double> times{0.0, 0.5};
  odelia::ode::OdeControl control(1e-6, 1e-6, 1.0, 0.0, 1e-12, 10.0, 1e-3);

  {
    CallbackSystem sys(lorenz_rhs, y0, 0.0, lorenz_jac, CallbackSystem::valid_type(), true);
    odelia::ode::Solver<CallbackSystem> s(sys, control, Method::rodas);
    s.advance_adaptive(times);
    const size_t accepted = s.times().size() - 1;
    const CallbackSystem& after = s.get_system_ref();
    std::printf("       (autonomous, analytic J: %zu accepted, %zu rejected, %zu rhs, %zu jac)\n",
                accepted, s.get_n_rejections(), after.n_rhs, after.n_jac);
    const size_t attempts = accepted + s.get_n_rejections();
    check(after.n_rhs == 1 + 6 * attempts,
          "an autonomous system with an analytic Jacobian costs 6 evaluations per attempt");
    check(after.n_jac == accepted,
          "and one Jacobian per accepted step, however many attempts it took");
  }
  {
    CallbackSystem sys(lorenz_rhs, y0, 0.0, lorenz_jac, CallbackSystem::valid_type(), false);
    odelia::ode::Solver<CallbackSystem> s(sys, control, Method::rodas);
    s.advance_adaptive(times);
    const size_t accepted = s.times().size() - 1;
    const size_t attempts = accepted + s.get_n_rejections();
    const CallbackSystem& after = s.get_system_ref();
    check(after.n_rhs == 1 + 6 * attempts + accepted,
          "a system not declared autonomous pays one more per accepted step for df/dt");
    check(after.n_jac == accepted, "and still one Jacobian per accepted step");
  }
  {
    CallbackSystem sys(lorenz_rhs, y0, 0.0, CallbackSystem::jac_type(), CallbackSystem::valid_type(), true);
    odelia::ode::Solver<CallbackSystem> s(sys, control, Method::rodas);
    s.advance_adaptive(times);
    const size_t accepted = s.times().size() - 1;
    const size_t attempts = accepted + s.get_n_rejections();
    const CallbackSystem& after = s.get_system_ref();
    check(after.n_rhs == 1 + 6 * attempts + 3 * accepted,
          "a finite-difference Jacobian adds n evaluations per accepted step, none per retry");
  }
  {
    // Force rejections with a first step far too large for the tolerance: the
    // Jacobian at the step start must not be re-formed for the retry.
    odelia::ode::OdeControl greedy(1e-8, 1e-8, 1.0, 0.0, 1e-12, 10.0, 1.0);
    CallbackSystem sys(lorenz_rhs, y0, 0.0, lorenz_jac, CallbackSystem::valid_type(), true);
    odelia::ode::Solver<CallbackSystem> s(sys, greedy, Method::rodas);
    s.advance_adaptive(times);
    const size_t accepted = s.times().size() - 1;
    const CallbackSystem& after = s.get_system_ref();
    std::printf("       (greedy start: %zu accepted, %zu rejected, %zu jac)\n",
                accepted, s.get_n_rejections(), after.n_jac);
    check(s.get_n_rejections() > 0, "a greedy first step is rejected (test is not vacuous)");
    check(after.n_jac == accepted, "and the Jacobian is reused across the retry");
  }
  {
    // The explicit stepper's budget is untouched by any of this: six per step.
    CallbackSystem sys(lorenz_rhs, y0, 0.0);
    odelia::ode::Solver<CallbackSystem> s(sys, control, Method::rkck);
    s.advance_adaptive(times);
    const size_t accepted = s.times().size() - 1;
    check(s.get_system_ref().n_rhs == 1 + 6 * accepted + s.get_n_rejections() * 6,
          "RKCK costs six evaluations per attempt");
  }
}

// A callable refusing a state by throwing DomainError costs the step, not the
// solve -- the #55 bargain, reached through a std::function.
void test_callback_domain_error_is_a_rejection() {
  auto logistic = [](double, const State& y, State& dydt) {
    if (y[0] < 0.0 || y[0] > 1.0) {
      odelia::util::stop_domain("y = " + odelia::util::format_double(y[0]) +
                                " is outside [0, 1]");
    }
    dydt[0] = 50.0 * y[0] * (1.0 - y[0]);
  };
  CallbackSystem sys(logistic, State{0.5}, 0.0);
  odelia::ode::Solver<CallbackSystem> s(sys, loose_control());

  bool threw = false;
  std::string msg;
  try {
    s.advance_adaptive(std::vector<double>{0.0, 1.0});
  } catch (const std::runtime_error& e) {
    threw = true;
    msg = e.what();
  }
  check(!threw, "a DomainError from a callback is a rejection, not a fatal");
  if (threw) {
    std::printf("       (raised: %s)\n", msg.c_str());
    return;
  }
  check(s.get_n_rejections() > 0, "the callback actually refused a step (test is not vacuous)");
  const double y = s.state()[0];
  check(y >= 0.0 && y <= 1.0, "and the solve finishes inside the domain");

  // The validity predicate is the other route, and a predicate is consulted on
  // the state vector rather than inside the callable.
  int refusals = 0;
  auto valid = [&refusals](double, const State& y) {
    const bool ok = y[0] >= 0.0 && y[0] <= 1.0;
    if (!ok) ++refusals;
    return ok;
  };
  auto plain = [](double, const State& y, State& dydt) {
    dydt[0] = 50.0 * y[0] * (1.0 - y[0]);
  };
  CallbackSystem checked(plain, State{0.5}, 0.0, CallbackSystem::jac_type(), valid);
  odelia::ode::Solver<CallbackSystem> s2(checked, loose_control());
  s2.advance_adaptive(std::vector<double>{0.0, 1.0});
  check(refusals > 0 && s2.state()[0] >= 0.0 && s2.state()[0] <= 1.0,
        "a validity predicate callback is enforced on the committed state");

  // Anything else a callable throws is a bug and ends the solve with its own
  // message.
  auto buggy = [](double, const State&, State&) {
    odelia::util::stop("bug-shaped failure, must not be absorbed");
  };
  CallbackSystem bug(buggy, State{0.5}, 0.0);
  threw = false;
  try {
    // Constructing the solver already evaluates the right-hand side once, to
    // seed itself, so the throw can come from here.
    odelia::ode::Solver<CallbackSystem> s3(bug, loose_control());
    s3.advance_adaptive(std::vector<double>{0.0, 1.0});
  } catch (const std::runtime_error& e) {
    threw = true;
    msg = e.what();
  }
  check(threw && msg.find("bug-shaped") != std::string::npos,
        "a util::stop() from a callback still ends the solve with its message");
}

// A consumer that drives the solver a step at a time and changes the shape of
// its state between steps (#62): the state can be re-seeded at a new size, the
// step size read back and set, a step bounded by a time it must not pass, and a
// stale Jacobian never survives a resize.
void test_callback_single_step_driving() {
  using odelia::ode::Method;
  // dy_i/dt = -y_i: linear, so the controller never rejects and a step of any
  // size lands where it was asked to.
  auto decay = [](double, const State& y, State& dydt) {
    for (size_t i = 0; i < y.size(); ++i) dydt[i] = -y[i];
  };
  auto decay_jac = [](double, const State& y, const State&, State& J) {
    const size_t n = y.size();
    for (size_t i = 0; i < n; ++i) J[i * n + i] = -1.0;
  };
  odelia::ode::OdeControl control(1e-6, 1e-6, 1.0, 0.0, 1e-12, 10.0, 1e-2);
  CallbackSystem sys(decay, State{1.0, 2.0, 3.0}, 0.0, decay_jac, CallbackSystem::valid_type(), true);
  odelia::ode::Solver<CallbackSystem> s(sys, control, Method::rodas);

  s.step(std::numeric_limits<double>::infinity());
  check(s.times().size() == 2 && s.times()[1] > 0.0, "one unbounded step advances the clock");
  check(s.state().size() == 3, "with three unknowns");

  // Shrink to two unknowns and carry on from a chosen state and time.
  const double t1 = s.time();
  s.get_system_ref().resize(2);
  s.set_state(State{0.5, 0.25}, t1);
  check(s.state().size() == 2 && s.time() == t1, "set_state re-seeds at a new size");
  check(s.get_step_size() == control.step_size_initial,
        "and the step size is back at the control's initial value");
  s.set_step_size(0.125);
  check(s.get_step_size() == 0.125, "set_step_size overrides it");
  s.step(std::numeric_limits<double>::infinity());
  check(std::abs((s.time() - t1) - 0.125) < 1e-15, "and the next step is exactly that size");
  check(std::abs(s.state()[0] - 0.5 * std::exp(-0.125)) < 1e-7,
        "with the right answer for the two-unknown system");

  // Grow to four, bounded steps: a step must stop at time_max, not pass it, and
  // the one that reaches it lands there exactly. (A step of 0.5 is rejected on
  // accuracy even for a linear system -- the method is fourth order, not exact
  // -- so this takes as many steps as the controller wants.)
  const double t2 = s.time();
  s.get_system_ref().resize(4);
  s.set_state(State{1.0, 1.0, 1.0, 1.0}, t2);
  s.set_step_size(10.0);
  int bounded_steps = 0;
  bool passed = false;
  while (s.time() < t2 + 0.5 && bounded_steps < 1000) {
    s.step(t2 + 0.5);
    passed = passed || s.time() > t2 + 0.5;
    ++bounded_steps;
  }
  check(!passed, "no bounded step passes time_max");
  check(s.time() == t2 + 0.5, "and the last one lands exactly on it");
  check(s.state().size() == 4, "with four unknowns");

  bool threw = false;
  std::string msg;
  try {
    s.step(t2 + 0.5);
  } catch (const std::runtime_error& e) {
    threw = true;
    msg = e.what();
  }
  check(threw && msg.find("already at time_max") != std::string::npos,
        "stepping from time_max itself is refused, not stepped by zero");

  // The last right-hand-side call of an accepted step is at the accepted state:
  // what a consumer keeping extra results from its last evaluation relies on.
  State last_y;
  double last_t = -1.0;
  auto recording = [&last_y, &last_t](double t, const State& y, State& dydt) {
    last_y = y;
    last_t = t;
    dydt[0] = -y[0];
  };
  CallbackSystem rec(recording, State{1.0}, 0.0, CallbackSystem::jac_type(), CallbackSystem::valid_type(), true);
  odelia::ode::Solver<CallbackSystem> sr(rec, control, Method::rodas);
  bool contract = true;
  for (int i = 0; i < 5; ++i) {
    sr.step(std::numeric_limits<double>::infinity());
    contract = contract && last_y == sr.state() && last_t == sr.time();
  }
  check(contract, "after each accepted RODAS step the last evaluation was at the accepted state");
  CallbackSystem rec2(recording, State{1.0}, 0.0);
  odelia::ode::Solver<CallbackSystem> sr2(rec2, control, Method::rkck);
  contract = true;
  for (int i = 0; i < 5; ++i) {
    sr2.step(std::numeric_limits<double>::infinity());
    contract = contract && last_y == sr2.state() && last_t == sr2.time();
  }
  check(contract, "and the same under RKCK");
}

// W = I/(h*gamma) - J can be singular at an unlucky h. That is a reason to try
// another step size, not to end the solve.
void test_singular_w_is_a_rejection() {
  using odelia::ode::Method;
  const double h0 = 1e-2, gamma = 0.25;
  int calls = 0;
  auto decay = [](double, const State& y, State& dydt) { dydt[0] = -y[0]; };
  // Exact Jacobian except on the first formation, which is rigged to make W
  // exactly zero at the control's initial step.
  auto jac = [&calls, h0, gamma](double, const State&, const State&, State& J) {
    ++calls;
    J[0] = (calls == 1) ? 1.0 / (h0 * gamma) : -1.0;
  };
  odelia::ode::OdeControl control(1e-6, 1e-6, 1.0, 0.0, 1e-12, 10.0, h0);
  CallbackSystem sys(decay, State{1.0}, 0.0, jac, CallbackSystem::valid_type(), true);
  odelia::ode::Solver<CallbackSystem> s(sys, control, Method::rodas);

  bool threw = false;
  std::string msg;
  try {
    s.advance_adaptive(std::vector<double>{0.0, 0.1});
  } catch (const std::runtime_error& e) {
    threw = true;
    msg = e.what();
  }
  check(!threw, "a singular W is a rejection, not a fatal");
  if (threw) {
    std::printf("       (raised: %s)\n", msg.c_str());
    return;
  }
  check(s.get_n_rejections() >= 1, "and it cost a retry");
  check(std::abs(s.state()[0] - std::exp(-0.1)) < 1e-5, "after which the solve is right");
}


// --- Dense output (#24) ------------------------------------------------------
//
// Landing a step on every requested output time costs steps the controller did
// not want: Lorenz to t = 100 at 1e-6 takes 27k evaluations with two output
// rows and 60k with ten thousand. The interpolant reads the state inside the
// last accepted step from its endpoints and their derivatives, which are
// already in hand, so a dense collect costs what the integration costs.

void test_dense_output() {
  using odelia::ode::Method;
  // dy/dt = 3 t^2, y = t^3: a cubic, which the Hermite interpolant reproduces
  // exactly, whatever the steps.
  auto cubic = [](double t, const State&, State& dydt) { dydt[0] = 3.0 * t * t; };
  odelia::ode::OdeControl control(1e-6, 1e-6, 1.0, 0.0, 1e-12, 0.3, 0.1);
  CallbackSystem sys(cubic, State{0.0}, 0.0);
  odelia::ode::Solver<CallbackSystem> s(sys, control, Method::rkck);
  std::vector<double> times;
  for (int i = 0; i <= 100; ++i) times.push_back(0.01 * i);
  const auto rows = s.advance_collect(times, true);
  double worst = 0.0;
  for (size_t k = 0; k < times.size(); ++k) {
    worst = std::max(worst, std::abs(rows[k][0] - times[k] * times[k] * times[k]));
  }
  check(rows.size() == times.size(), "a dense collect returns one row per time");
  check(worst < 1e-12, "and reproduces a cubic exactly at every requested time");
  check(s.time() == times.back(), "landing exactly on the last time");
  const size_t steps = s.times().size() - 1;
  check(steps < times.size() - 1, "with fewer steps than output times (test is not vacuous)");
  std::printf("       (%zu steps for %zu output rows, worst error %.1e)\n",
              steps, times.size(), worst);

  // Lorenz: the dense rows agree with landed rows to within the tolerance's
  // reach, for many fewer evaluations.
  const State y0{1.0, 1.0, 1.0};
  std::vector<double> many;
  for (int i = 0; i <= 2000; ++i) many.push_back(0.001 * i);
  odelia::ode::OdeControl tight(1e-8, 1e-8, 1.0, 0.0, 1e-12, 10.0, 1e-6);
  CallbackSystem a(lorenz_rhs, y0, 0.0, lorenz_jac, CallbackSystem::valid_type(), true);
  odelia::ode::Solver<CallbackSystem> sa(a, tight, Method::rkck);
  const auto dense = sa.advance_collect(many, true);
  CallbackSystem b(lorenz_rhs, y0, 0.0, lorenz_jac, CallbackSystem::valid_type(), true);
  odelia::ode::Solver<CallbackSystem> sb(b, tight, Method::rkck);
  const auto landed = sb.advance_collect(many, false);
  double worst_l = 0.0;
  for (size_t k = 0; k < many.size(); ++k) {
    worst_l = std::max(worst_l, max_abs_diff(dense[k], landed[k]));
  }
  check(worst_l < 1e-4, "dense and landed Lorenz rows agree to the interpolant's order");
  check(sa.get_system_ref().n_rhs < sb.get_system_ref().n_rhs / 2,
        "and the dense collect made under half the evaluations");
  std::printf("       (dense %zu evaluations, landed %zu, worst difference %.1e)\n",
              sa.get_system_ref().n_rhs, sb.get_system_ref().n_rhs, worst_l);

  // The same under RODAS, whose end derivative is now carried forward.
  CallbackSystem c(lorenz_rhs, y0, 0.0, lorenz_jac, CallbackSystem::valid_type(), true);
  odelia::ode::Solver<CallbackSystem> sc(c, tight, Method::rodas);
  const auto dense_r = sc.advance_collect(many, true);
  double worst_r = 0.0;
  for (size_t k = 0; k < many.size(); ++k) {
    worst_r = std::max(worst_r, max_abs_diff(dense_r[k], landed[k]));
  }
  check(worst_r < 1e-4, "and RODAS dense rows agree with landed RKCK rows");

  // Refusals: no accepted step, a time outside the last step, bad times.
  CallbackSystem d(cubic, State{0.0}, 0.0);
  odelia::ode::Solver<CallbackSystem> sd(d, control, Method::rkck);
  bool threw = false;
  try { State out; sd.get_internal().interpolate(0.0, out); } catch (const std::runtime_error&) { threw = true; }
  check(threw, "interpolating before any accepted step is refused");
  sd.step(0.05);
  threw = false;
  try { State out; sd.get_internal().interpolate(0.06, out); } catch (const std::runtime_error&) { threw = true; }
  check(threw, "and so is a time outside the last step");
  threw = false;
  try { sd.advance_collect(std::vector<double>{sd.time(), 1.0, 0.5}, true); } catch (const std::runtime_error&) { threw = true; }
  check(threw, "and times that are not increasing");
}


// --- Dormand-Prince 5(4) with order-4 dense output -----------------------------

void test_dopri() {
  using odelia::ode::Method;
  const State y0{1.0, 1.0, 1.0};
  const std::vector<double> times{0.0, 1.0, 2.0};

  // Agrees with the other steppers.
  CallbackSystem a(lorenz_rhs, y0, 0.0);
  odelia::ode::Solver<CallbackSystem> sa(a, tight_control(), Method::rkck);
  sa.advance_adaptive(times);
  CallbackSystem b(lorenz_rhs, y0, 0.0);
  odelia::ode::Solver<CallbackSystem> sb(b, tight_control(), Method::dopri);
  sb.advance_adaptive(times);
  check(max_abs_diff(sa.state(), sb.state()) < 1e-5, "Dormand-Prince agrees with RKCK on Lorenz");

  // Six evaluations per attempt, first-same-as-last.
  odelia::ode::OdeControl control(1e-6, 1e-6, 1.0, 0.0, 1e-12, 10.0, 1e-3);
  CallbackSystem c(lorenz_rhs, y0, 0.0);
  odelia::ode::Solver<CallbackSystem> sc(c, control, Method::dopri);
  sc.advance_adaptive(std::vector<double>{0.0, 0.5});
  const size_t attempts = sc.times().size() - 1 + sc.get_n_rejections();
  check(sc.get_system_ref().n_rhs == 1 + 6 * attempts, "and costs six evaluations per attempt");

  // The dense output is exact for a quartic: y' = 4 t^3, y = t^4. A wrong
  // dense coefficient shows up here at O(1e-3).
  auto quartic = [](double t, const State&, State& dydt) { dydt[0] = 4.0 * t * t * t; };
  odelia::ode::OdeControl coarse(1e-6, 1e-6, 1.0, 0.0, 1e-12, 0.4, 0.25);
  CallbackSystem q(quartic, State{0.0}, 0.0);
  odelia::ode::Solver<CallbackSystem> sq(q, coarse, Method::dopri);
  std::vector<double> grid;
  for (int i = 0; i <= 200; ++i) grid.push_back(0.01 * i);
  const auto rows = sq.advance_collect(grid, true);
  double worst = 0.0;
  for (size_t k = 0; k < grid.size(); ++k) {
    const double t = grid[k];
    worst = std::max(worst, std::abs(rows[k][0] - t * t * t * t));
  }
  check(worst < 1e-12, "its dense output reproduces a quartic exactly");
  check(sq.times().size() - 1 < 20, "from a handful of steps (test is not vacuous)");
  std::printf("       (%zu steps for %zu rows, worst error %.1e)\n", sq.times().size() - 1, grid.size(), worst);

  // On Lorenz the dense rows are within the tolerance's reach of the landed
  // rows, where cubic Hermite was 40-100x out.
  std::vector<double> many;
  for (int i = 0; i <= 2000; ++i) many.push_back(0.001 * i);
  odelia::ode::OdeControl tol8(1e-8, 1e-8, 1.0, 0.0, 1e-12, 10.0, 1e-6);
  CallbackSystem d(lorenz_rhs, y0, 0.0);
  odelia::ode::Solver<CallbackSystem> sd(d, tol8, Method::dopri);
  const auto dense = sd.advance_collect(many, true);
  CallbackSystem e(lorenz_rhs, y0, 0.0);
  odelia::ode::Solver<CallbackSystem> se(e, tol8, Method::dopri);
  const auto landed = se.advance_collect(many, false);
  double worst_l = 0.0;
  for (size_t k = 0; k < many.size(); ++k) worst_l = std::max(worst_l, max_abs_diff(dense[k], landed[k]));
  check(worst_l < 1e-5, "dense Lorenz rows are within 1e-5 of landed rows at 1e-8");
  check(sd.get_system_ref().n_rhs < se.get_system_ref().n_rhs / 2, "for under half the evaluations");
  std::printf("       (dense %zu evaluations, landed %zu, worst difference %.1e)\n",
              sd.get_system_ref().n_rhs, se.get_system_ref().n_rhs, worst_l);
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
  test_pinned_step_domain_error_is_a_rejection();
  test_pinned_step_predicate_is_enforced();
  test_pinned_step_non_domain_throw_is_not_absorbed();
  test_pinned_step_unreachable_domain_fails_with_a_reason();
  test_pinned_step_unchanged_when_nothing_objects();
  test_jacobian_hook_is_detected();
  test_callback_rodas_matches_compiled();
  test_callback_call_budget();
  test_callback_domain_error_is_a_rejection();
  test_callback_single_step_driving();
  test_singular_w_is_a_rejection();
  test_dense_output();
  test_dopri();
  if (failures == 0) {
    std::printf("all checks passed\n");
    return 0;
  }
  std::printf("%d failure(s)\n", failures);
  return 1;
}
