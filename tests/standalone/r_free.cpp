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
// includes and only receives by accident from R's headers. A retired spline
// header came to use `assert` without <cassert>, which is how.
//
// The R interface headers -- solver_interface.hpp, rcpp_interface_helpers.hpp
// -- are deliberately NOT listed here. They are meant to depend on Rcpp.
//
// Note that this links src/Tape.cpp: the XAD tape runtime lives in exactly one
// object file, by design (see ARCHITECTURE.md), so a standalone consumer that
// instantiates Solver has to compile it too. That is a linking requirement, not
// an R one.

#include <odelia/ode_util.hpp>
#include <odelia/interpolator.hpp>
#include <odelia/drivers.hpp>
#include <odelia/ode_control.hpp>
#include <odelia/ode_interface.hpp>
#include <odelia/ode_step_rkck.hpp>
#include <odelia/ode_step_rodas.hpp>
#include <odelia/ode_step_dopri.hpp>
#include <odelia/ode_solver_internal.hpp>
#include <odelia/ode_solver.hpp>
#include <odelia/sweep.hpp>
#include <odelia/implicit_node.hpp>
#include <odelia/tangent.hpp>
#include <odelia/ode_steady_state.hpp>
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
  odelia::interpolator::hermite_interpolator<double> in;
  // x^2 with its own slope at every knot, which a cubic reproduces exactly.
  in.init({0.0, 1.0, 2.0, 3.0}, {0.0, 1.0, 4.0, 9.0}, {0.0, 2.0, 4.0, 6.0});
  check(std::abs(in.eval(2.0) - 4.0) < 1e-12, "interpolator hits its knots");
  check(std::abs(in.eval(1.5) - 2.25) < 1e-12, "and the quadratic between them");
  check(std::abs(in.slope(1.5) - 3.0) < 1e-12, "with the slope of the same curve");

  bool threw = false;
  try {
    odelia::interpolator::hermite_interpolator<double> one_knot;
    one_knot.init({0.0}, {0.0}, {0.0});
  } catch (const std::runtime_error &) {
    threw = true;
  }
  check(threw, "interpolator rejects a knot set with no span");
}

// A driver given as values alone: natural by default, monotone on request. An
// intermittent non-negative series is the case the choice exists for -- a single
// wet day between dry ones pulls a natural spline below zero beside it.
void test_driver_slopes() {
  const std::vector<double> x{0, 1, 2, 3, 4, 5, 6};
  const std::vector<double> y{0, 0, 0, 8, 0, 0, 0};
  odelia::drivers::Drivers d;
  d.set_variable("natural", x, y);
  d.set_variable("monotone", x, y, odelia::drivers::Slopes::monotone);
  double lo_nat = 0.0, lo_mono = 0.0, hi_mono = 0.0;
  for (double u = 0.0; u <= 6.0; u += 0.01) {
    lo_nat = std::min(lo_nat, d.evaluate("natural", u));
    lo_mono = std::min(lo_mono, d.evaluate("monotone", u));
    hi_mono = std::max(hi_mono, d.evaluate("monotone", u));
  }
  check(lo_nat < 0.0, "a natural driver dips below an intermittent series");
  check(lo_mono >= 0.0 && hi_mono <= 8.0,
        "a monotone driver stays inside the values bracketing each span");
  check(d.evaluate("monotone", 3.0) == 8.0, "and still hits its knots");
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
// template parameter so that the silent copy genuinely lacks the method and
// ChecksState resolves to false for it -- an `if constexpr` inside one
// struct would still leave the member there for the concept to find.
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
  using domain = std::vector<double>;
  check(odelia::ode::ChecksState<LogisticChecked, domain>,
        "ChecksState finds a declared ode_state_valid");
  check(!odelia::ode::ChecksState<Logistic<OnLeave::nothing>, domain>,
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
  std::printf("       (%d refusal(s); unguarded copy finished at y = %g)\n",
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


// A supplied row set is ONE statement, whatever the row count, and the rows it
// carries are the ones it was handed. The count is the guard: written as a sum of
// `out += d * (x - to_passive(x))` the same call is one recorded assignment per
// row, which is how a submodel's whole arithmetic reaches a consumer's tape.
void test_supplied_rows_cost_one_statement() {
  using A = odelia::ode::active_scalar<double>;
  using Tape = odelia::ode::adjoint_tape<double>;

  for (int n : {1, 5, 31}) {
    Tape tape;
    std::vector<A> x(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
      x[static_cast<std::size_t>(i)] = 1.0 + 0.25 * double(i);
    }
    // Never registered, so it holds no slot: its row must be dropped rather than
    // pushed, and the sweep must survive it.
    const A unregistered = 9.0;
    tape.registerInputs(x.begin(), x.end());
    tape.newRecording();

    std::vector<odelia::input_and_derivative<A>> against;
    for (int i = 0; i < n; ++i) {
      against.push_back({x[static_cast<std::size_t>(i)], 1.0 / (double(i) + 2.0)});
    }
    against.push_back({unregistered, 4.0});
    against.push_back({x[0], 0.0});

    const std::size_t s0 = tape.getNumStatements();
    A out;
    const odelia::record_report report =
        odelia::record_with_derivatives<A>(7.5, against, out);
    const std::size_t statements = tape.getNumStatements() - s0;

    tape.registerOutput(out);
    xad::derivative(out) = 1.0;
    tape.computeAdjoints();

    double worst = std::fabs(xad::value(out) - 7.5);
    for (int i = 0; i < n; ++i) {
      worst = std::fmax(worst, std::fabs(xad::derivative(x[static_cast<std::size_t>(i)]) -
                                         1.0 / (double(i) + 2.0)));
    }
    const std::string at = " (" + std::to_string(n) + " rows)";
    check(report.whole, "every row is recorded" + at);
    check(statements == 1, "one statement" + at);
    check(worst < 1e-15, "the value and every row are the ones supplied" + at);
  }
}

// The same call at a direction, which has no tape to hold a statement: the rows
// are the arithmetic there, and dropping them would be silent.
void test_supplied_rows_carry_a_direction() {
  using T = odelia::ode::tangent_scalar<double>;
  T x = 2.0;
  odelia::ode::seed_direction(x, 1.0);
  std::vector<odelia::input_and_derivative<T>> against{{x, 0.25}};
  T out;
  const odelia::record_report report =
      odelia::record_with_derivatives<T>(7.5, against, out);
  check(report.whole, "a direction records its rows");
  check(std::fabs(odelia::util::to_passive(out) - 7.5) < 1e-15,
        "the value is untouched at a direction");
  check(std::fabs(odelia::ode::derivative_along(out) - 0.25) < 1e-15,
        "and the direction carries the supplied row");
}


// A residual taken off the caller's tape gives the same rows as one left on it.
//
// F(p) = p^2 - x*y has the root p* = sqrt(x*y), so dp*/dx = y/(2p*) and
// dp*/dy = x/(2p*) in closed form -- which is the referee here, rather than the
// two routes agreeing with each other.
void test_a_preaccumulated_residual_keeps_its_rows() {
  using A = odelia::ode::active_scalar<double>;
  using Tape = odelia::ode::adjoint_tape<double>;
  const double x0 = 2.0, y0 = 8.0;
  const double root = 4.0;          // sqrt(2*8)
  const double dFdp = 2.0 * root;   // 8
  const double want_dx = y0 / dFdp; // 1.0
  const double want_dy = x0 / dFdp; // 0.25

  double got_dx[2], got_dy[2];
  std::size_t statements[2];


  for (int arm = 0; arm < 2; ++arm) {
    Tape tape;
    A x = x0, y = y0;
    tape.registerInput(x);
    tape.registerInput(y);
    tape.newRecording();
    auto residual = [&](const A& p) -> A { return p * p - x * y; };
    const std::size_t s0 = tape.getNumStatements();
    A p_star = (arm == 0)
                   ? odelia::implicit_value<A>(root, dFdp, residual)
                   : odelia::implicit_value<A>(root, dFdp, residual, x, y);
    statements[arm] = tape.getNumStatements() - s0;
    tape.registerOutput(p_star);
    xad::derivative(p_star) = 1.0;
    tape.computeAdjoints();
    got_dx[arm] = xad::derivative(x);
    got_dy[arm] = xad::derivative(y);
    check(std::fabs(xad::value(p_star) - root) < 1e-14,
          arm == 0 ? "the value is the root (on the tape)"
                   : "the value is the root (preaccumulated)");
    if (arm == 1) {
      check(odelia::ode::count_active_slots<A>(x, y) == 2,
            "the walk had both inputs to reach");
    }
  }

  check(std::fabs(got_dx[0] - want_dx) < 1e-12 &&
            std::fabs(got_dy[0] - want_dy) < 1e-12,
        "the recorded residual gives the theorem's rows");
  check(std::fabs(got_dx[1] - want_dx) < 1e-12 &&
            std::fabs(got_dy[1] - want_dy) < 1e-12,
        "and so does the preaccumulated one");
  check(statements[1] == 1, "which costs one statement");
  check(statements[1] < statements[0],
        "against the whole residual left on the tape");
  std::printf("       (on the tape %zu statements, preaccumulated %zu)\n",
              statements[0], statements[1]);
}

// The rows accumulate, so a second solve against the same inputs must not add to
// the first. Every number stays finite when it does, which is what makes it worth
// a check of its own.
void test_two_preaccumulated_solves_do_not_add_up() {
  using A = odelia::ode::active_scalar<double>;
  using Tape = odelia::ode::adjoint_tape<double>;
  Tape tape;
  A x = 2.0, y = 8.0;
  tape.registerInput(x);
  tape.registerInput(y);
  tape.newRecording();
  auto residual = [&](const A& p) -> A { return p * p - x * y; };
  const A first = odelia::implicit_value<A>(4.0, 8.0, residual, x, y);
  const A second = odelia::implicit_value<A>(4.0, 8.0, residual, x, y);
  A sum = first + second;
  tape.registerOutput(sum);
  xad::derivative(sum) = 1.0;
  tape.computeAdjoints();
  // Two identical solves, so each row is twice one solve's and no more.
  check(std::fabs(xad::derivative(x) - 2.0 * 1.0) < 1e-12 &&
            std::fabs(xad::derivative(y) - 2.0 * 0.25) < 1e-12,
        "two solves carry two rows, not three");
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
// --- A recorded run, swept and replayed -------------------------------------
//
// y' = r y (1 - y) on [0, 1], with the domain declared. Carries what a sweep
// needs (rebind_from, ad_parameters, for_each_active, set_recorded_state), so
// solve_adjoint can be refereed here against a central difference of the same
// pinned run.
namespace {

template <typename T = double>
struct Grow {
  using value_type = T;
  template <typename> friend struct Grow;
  T r, y, dydt, y_init;
  double time = 0.0, t0 = 0.0;
  explicit Grow(T r_ = T(1.0), double y0 = 0.5)
      : r(r_), y(y0), dydt(0.0), y_init(y0) { compute_rates(); }
  template <class S> Grow<S> rebind_from() const {
    Grow<S> g(S(xad::value(r)), xad::value(y_init));
    std::vector<S> st{S(xad::value(y))};
    g.set_ode_state(st.begin(), time);
    return g;
  }
  size_t ode_size() const { return 1; }
  double ode_time() const { return time; }
  template <class It> It set_ode_state(It it, double t) {
    y = *it++; time = t; compute_rates(); return it;
  }
  void set_recorded_state(const std::vector<T>& s, double t) {
    set_ode_state(s.begin(), t);
  }
  void compute_rates() { dydt = r * y * (T(1.0) - y); }
  template <class It> It ode_state(It it) const { *it++ = y; return it; }
  template <class It> It ode_rates(It it) const { *it++ = dydt; return it; }
  template <class It> It ode_initial_state(It it) const { *it++ = y_init; return it; }
  template <class It> It set_initial_state(It it, double t0_ = 0.0) {
    t0 = t0_; y_init = *it++; return it;
  }
  void reset() { y = y_init; time = t0; compute_rates(); }
  std::vector<T*> ad_parameters() { return {&r}; }
  template <class F> void for_each_active(F&& f) { f(r); f(y); f(dydt); f(y_init); }
  bool ode_state_valid(const std::vector<T>& s) const {
    return xad::value(s[0]) >= 0.0 && xad::value(s[0]) <= 1.0;
  }
};

double grow_pinned_final(double r, const std::vector<double>& times,
                         const odelia::ode::OdeControl& c) {
  Grow<double> g(r, 0.5);
  odelia::ode::Solver<Grow<double>> s(g, c);
  s.advance_fixed(times);
  return s.state()[0];
}

} // namespace

// Pinned finely enough that nothing is refused, the sweep of a pinned run is the
// derivative of that run: refereed against a central difference, not against
// itself.
void test_sweep_of_a_pinned_run_matches_a_difference() {
  const std::vector<double> times{0.0, 0.01, 0.02, 0.03, 0.04, 0.05};
  const double r = 50.0;
  Grow<double> g(r, 0.5);
  odelia::ode::Solver<Grow<double>> s(g, loose_control());
  s.set_keep_states(true);
  s.reset();
  s.advance_fixed(times);
  check(s.get_n_rejections() == 0, "no pinned interval was refused");
  check(s.recording().size() == times.size(), "one row per pinned time");

  odelia::ode::adjoint_rows lambda = odelia::ode::adjoint_rows::one_row({1.0});
  odelia::ode::adjoint_rows dp(1, 1);
  s.solve_adjoint(lambda, dp);
  const double h = 1e-6;
  const double fd = (grow_pinned_final(r + h, times, loose_control()) -
                     grow_pinned_final(r - h, times, loose_control())) / (2 * h);
  check(std::fabs(dp[0][0] - fd) < 1e-7 * std::fabs(fd),
        "d y(T)/dr from the sweep is the central difference of the run");
  std::printf("       (sweep %.10g, difference %.10g)\n", dp[0][0], fd);
}

// A pinned interval the domain refused is crossed in several sub-steps and
// recorded as one row. Swept as one step of the interval it gave a derivative
// forty orders of magnitude off, finitely; now the row is marked and refused.
void test_subdivided_pinned_row_is_refused() {
  const std::vector<double> times{0.0, 0.5, 1.0};
  Grow<double> g(50.0, 0.5);
  odelia::ode::Solver<Grow<double>> s(g, loose_control());
  s.set_keep_states(true);
  s.reset();
  s.advance_fixed(times);
  check(s.get_n_rejections() > 0, "at least one pinned interval was refused");
  const auto rec = s.recording();
  check(rec.size() == times.size(), "the record still holds one row per time");
  bool marked = false;
  for (const auto& row : rec) marked = marked || row.subdivided;
  check(marked, "and the row that was subdivided says so");

  odelia::ode::adjoint_rows lambda = odelia::ode::adjoint_rows::one_row({1.0});
  odelia::ode::adjoint_rows dp(1, 1);
  bool refused = false;
  try {
    s.solve_adjoint(lambda, dp);
  } catch (const std::runtime_error& e) {
    refused = std::string(e.what()).find("sub-steps") != std::string::npos;
  }
  check(refused, "the sweep refuses the subdivided row by name");

  Grow<double> again(50.0, 0.5);
  odelia::ode::Solver<Grow<double>> replay(again, loose_control());
  refused = false;
  try {
    replay.advance_recorded(rec);
  } catch (const std::runtime_error& e) {
    refused = std::string(e.what()).find("subdivided") != std::string::npos;
  }
  check(refused, "and so does a replay of the recording");
}

// The two batches a sweep takes are written differently -- one replaced, one
// accumulated -- so the same object for both loses whichever was written first.
void test_sweep_refuses_one_batch_for_both() {
  LorenzSystem<double> sys(10.0, 28.0, 8.0 / 3.0);
  odelia::ode::Solver<LorenzSystem<double>> s(sys, odelia::ode::OdeControl());
  s.set_keep_states(true);
  s.reset();
  s.advance_adaptive(std::vector<double>{0.0, 0.5});
  odelia::ode::adjoint_rows lambda =
      odelia::ode::adjoint_rows::one_row({1.0, 0.0, 0.0});
  bool refused = false;
  try {
    s.solve_adjoint(lambda, lambda);
  } catch (const std::runtime_error& e) {
    refused = std::string(e.what()).find("same batch") != std::string::npos;
  }
  check(refused, "solve_adjoint refuses the same batch as seed and accumulator");
}

// An input that already carries an adjoint when a residual is taken off the
// tape must get exactly that adjoint back and a row that does not include it.
// With held adjoints of zero the two are indistinguishable, which is why this
// is a case of its own.
void test_implicit_value_leaves_the_callers_adjoint_alone() {
  using A = odelia::ode::active_scalar<double>;
  using Tape = odelia::ode::adjoint_tape<double>;
  const double x0 = 2.0, y0 = 8.0, held = 9.0;
  const double root = 4.0, dFdp = 2.0 * root;
  Tape tape;
  A x = x0, y = y0;
  tape.registerInput(x);
  tape.registerInput(y);
  tape.newRecording();
  xad::derivative(x) = held;
  xad::derivative(y) = held;
  auto residual = [&](const A& p) -> A { return p * p - x * y; };
  A p_star = odelia::implicit_value<A>(root, dFdp, residual, x, y);
  check(xad::derivative(x) == held && xad::derivative(y) == held,
        "the inputs' adjoints are as the caller left them");
  tape.registerOutput(p_star);
  xad::derivative(p_star) = 1.0;
  tape.computeAdjoints();
  // The theorem's rows on top of what was held: dp*/dx = y/(2p*), dp*/dy = x/(2p*).
  check(std::fabs(xad::derivative(x) - (held + y0 / dFdp)) < 1e-12 &&
            std::fabs(xad::derivative(y) - (held + x0 / dFdp)) < 1e-12,
        "the rows add the theorem's derivative and nothing of the held adjoint");
}

// --- A run whose state widens -------------------------------------------------
//
// Cells of y' = r y (1 - y), one more inserted mid-run with a size set by r, so
// the newborn's initial condition is a channel of parameter derivative of its
// own, beside the dynamics'.
namespace {

template <typename T = double>
struct Cells {
  using value_type = T;
  template <typename> friend struct Cells;
  T r;
  std::vector<T> y, dydt;
  std::vector<double> y_init;
  double time = 0.0;
  explicit Cells(T r_ = T(1.0), std::vector<double> y0 = {0.5})
      : r(r_), y_init(std::move(y0)) { reset(); }
  template <class S> Cells<S> rebind_from() const {
    std::vector<double> v;
    for (const T& c : y) v.push_back(xad::value(c));
    Cells<S> out(S(xad::value(r)), y_init);
    out.y.assign(v.begin(), v.end());
    out.dydt.resize(v.size());
    out.time = time;
    out.compute_rates();
    return out;
  }
  size_t ode_size() const { return y.size(); }
  double ode_time() const { return time; }
  template <class It> It set_ode_state(It it, double t) {
    for (T& c : y) c = *it++;
    time = t; compute_rates(); return it;
  }
  void set_recorded_state(const std::vector<T>& s, double t) {
    y.resize(s.size()); dydt.resize(s.size());
    set_ode_state(s.begin(), t);
  }
  void compute_rates() {
    for (size_t i = 0; i < y.size(); ++i) dydt[i] = r * y[i] * (T(1.0) - y[i]);
  }
  template <class It> It ode_state(It it) const { for (const T& c : y) *it++ = c; return it; }
  template <class It> It ode_rates(It it) const { for (const T& c : dydt) *it++ = c; return it; }
  void reset() {
    y.assign(y_init.begin(), y_init.end()); dydt.resize(y.size());
    time = 0.0; compute_rates();
  }
  std::vector<T*> ad_parameters() { return {&r}; }
  template <class F> void for_each_active(F&& f) {
    f(r); for (T& c : y) f(c); for (T& c : dydt) f(c);
  }
  // The cells as they were, then a newborn of size r / 100.
  template <class It> void apply_insertion(double, It x, std::vector<T>& out) {
    const size_t n = y.size();
    y.resize(n + 1); dydt.resize(n + 1);
    for (size_t i = 0; i < n; ++i) y[i] = *x++;
    y[n] = T(0.01) * r;
    compute_rates();
    out.resize(n + 1);
    ode_state(out.begin());
  }
};

// Widen the solver's System at its current time and record the row.
void insert_cell(odelia::ode::Solver<Cells<double>>& s) {
  Cells<double>& sys = s.get_system_ref();
  std::vector<double> before(sys.ode_size());
  sys.ode_state(before.begin());
  std::vector<double> widened;
  sys.apply_insertion(s.time(), before.begin(), widened);
  s.set_state_from_system();
  s.push_insertion();
}

// The run: to 0.5, insert, to t_end. Pinned to `seg1`/`seg2` when given.
double cells_sum(double r, const std::vector<double>* seg1,
                 const std::vector<double>* seg2, double t_end = 1.0) {
  Cells<double> c(r, {0.5});
  odelia::ode::Solver<Cells<double>> s(c, odelia::ode::OdeControl());
  if (seg1) s.advance_fixed(*seg1); else s.advance_adaptive(std::vector<double>{0.0, 0.5});
  insert_cell(s);
  if (seg2) s.advance_fixed(*seg2); else s.advance_adaptive(std::vector<double>{0.5, t_end});
  double sum = 0.0;
  for (double v : s.state()) sum += v;
  return sum;
}

} // namespace

void test_sweep_across_an_insertion() {
  const double r = 2.0;
  Cells<double> c(r, {0.5});
  odelia::ode::Solver<Cells<double>> s(c, odelia::ode::OdeControl());
  s.set_keep_states(true);
  s.reset();
  s.advance_adaptive(std::vector<double>{0.0, 0.5});
  insert_cell(s);
  s.advance_adaptive(std::vector<double>{0.5, 1.0});
  const auto rec = s.recording();
  size_t insertions = 0;
  for (const auto& row : rec) insertions += row.insertion;
  check(insertions == 1, "the recording holds the insertion row");
  check(rec.front().state.size() == 1 && rec.back().state.size() == 2,
        "and widens from one cell to two");

  // The same steps, pinned, for the difference.
  std::vector<double> seg1, seg2;
  for (const auto& ins : s.schedule()) {
    if (ins.time <= 0.5) seg1.push_back(ins.time);
    if (ins.time >= 0.5) seg2.push_back(ins.time);
  }
  odelia::ode::adjoint_rows lambda = odelia::ode::adjoint_rows::one_row({1.0, 1.0});
  odelia::ode::adjoint_rows dp(1, 1);
  const size_t ranges = s.solve_adjoint(lambda, dp);
  check(ranges == 2, "two ranges, one each side of the insertion");
  check(lambda.width() == 1, "the adjoint comes back at the width the run started at");
  const double h = 1e-6;
  const double fd = (cells_sum(r + h, &seg1, &seg2) - cells_sum(r - h, &seg1, &seg2)) / (2 * h);
  check(std::fabs(dp[0][0] - fd) < 1e-7 * std::fabs(fd),
        "d(sum of cells)/dr across the insertion is the central difference");
  std::printf("       (sweep %.10g, difference %.10g)\n", dp[0][0], fd);

  // The newborn's size is r / 100, so its own row carries a derivative that is
  // not the dynamics'.
  odelia::ode::adjoint_rows newborn = odelia::ode::adjoint_rows::one_row({0.0, 1.0});
  odelia::ode::adjoint_rows dp_newborn(1, 1);
  s.solve_adjoint(newborn, dp_newborn);
  check(dp_newborn[0][0] != 0.0 && std::fabs(dp_newborn[0][0]) < std::fabs(dp[0][0]),
        "the newborn's initial condition contributes a parameter derivative");

  // After the sweep the solver stands where the run left it, at the run's width,
  // and steps on from there as a run that was never swept does.
  check(s.state().size() == 2, "the solver's state is at the run's width after the sweep");
  s.advance_adaptive(std::vector<double>{1.0, 1.5});
  Cells<double> c2(r, {0.5});
  odelia::ode::Solver<Cells<double>> unswept(c2, odelia::ode::OdeControl());
  unswept.advance_adaptive(std::vector<double>{0.0, 0.5});
  insert_cell(unswept);
  unswept.advance_adaptive(std::vector<double>{0.5, 1.0});
  unswept.advance_adaptive(std::vector<double>{1.0, 1.5});
  check(s.state() == unswept.state() && s.time() == unswept.time(),
        "and a step taken after the sweep is the step an unswept run takes");

  // A replay of the recording records the insertion row too, so its own
  // recording can be swept.
  Cells<double> again(r, {0.5});
  odelia::ode::Solver<Cells<double>> replay(again, odelia::ode::OdeControl());
  replay.set_keep_states(true);
  replay.reset();
  replay.advance_recorded(rec);
  const auto rec2 = replay.recording();
  check(rec2.size() == rec.size(), "a replay's recording has the run's rows");
  bool same = true;
  for (size_t k = 0; k < rec.size(); ++k) {
    same = same && rec2[k].insertion == rec[k].insertion &&
           rec2[k].state.size() == rec[k].state.size();
  }
  check(same, "with the insertion where the run had it");
}

// A scalar that already carries an adjoint is closed to the forward-AD Jacobian,
// and naming the Jacobian at it must not be a compile error: the implicit
// stepper refuses at run time, as it did before the tangent scalar was named.
void test_jacobian_is_closed_at_an_adjoint_scalar() {
  using A = odelia::ode::active_scalar<double>;
  using Sys = LorenzSystem<A>;
  check(odelia::ode::Jacobian<Sys>::value_is_adjoint,
        "the Jacobian sees an adjoint scalar");
  check(!odelia::ode::Jacobian<Sys>::ad_supported && !odelia::ode::Jacobian<Sys>::supported,
        "and offers no route to a Jacobian there");
  check(!odelia::ode::RodasStep<Sys>::supported, "so RODAS is closed to it");
  // The step-size controller reads doubles, so an adjoint-scalar run is pinned;
  // the adaptive path at such a scalar has never compiled.
  Sys sys(A(10.0), A(28.0), A(8.0 / 3.0));
  odelia::ode::Solver<Sys> s(sys, odelia::ode::OdeControl(), odelia::ode::Method::rodas);
  bool refused = false;
  try {
    s.advance_fixed(std::vector<double>{0.0, 0.1});
  } catch (const std::runtime_error& e) {
    refused = std::string(e.what()).find("not available") != std::string::npos;
  }
  check(refused, "and a RODAS step at an adjoint scalar is refused at run time");
  odelia::ode::Solver<Sys> rk(sys, odelia::ode::OdeControl());
  rk.advance_fixed(std::vector<double>{0.0, 0.1});
  check(rk.time() == 0.1, "while the explicit stepper steps at it as before");
}

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

// A compiled system with its own Jacobian hook and no rebind_from(): the hook alone
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
  using odelia::ode::Rebindable;
  using odelia::ode::tangent_scalar;
  using odelia::ode::Jacobian;
  using odelia::ode::RodasStep;

  check(has_jacobian<HookedLorenz>::value, "has_jacobian finds a declared ode_jacobian");
  check(!Rebindable<HookedLorenz, tangent_scalar<double>>, "on a system with no rebind_from()");
  check(Jacobian<HookedLorenz>::supported && !Jacobian<HookedLorenz>::ad_supported,
        "so the Jacobian is supported through the hook, not AD");
  check(RodasStep<HookedLorenz>::supported, "and RODAS is open to it");

  check(has_jacobian<CallbackSystem>::value, "CallbackSystem declares the hook");
  check(RodasStep<CallbackSystem>::supported, "and RODAS is open to it too");

  check(!has_jacobian<LorenzSystem<double>>::value &&
            Jacobian<LorenzSystem<double>>::ad_supported,
        "the compiled Lorenz keeps its AD route (no hook, has rebind_from)");
  check(!Jacobian<Logistic<OnLeave::nothing>>::supported,
        "a system with neither hook nor rebind_from is not supported");

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

  // The default floor keeps the perturbation above rounding: f = exp(y) + 1
  // has unit slope at y = 0, where a floor of a tight absolute tolerance
  // (1e-10) makes the perturbation 1e-16 and the difference vanishes.
  auto expo = [](double, const State& y, State& dydt) { dydt[0] = std::exp(y[0]) + 1.0; };
  CallbackSystem ex(expo, State{0.0}, 0.0);
  State y_zero{0.0}, f_zero(1), J1, J2;
  odelia::ode::derivs(ex, y_zero, f_zero, 0.0);
  odelia::ode::fd_jacobian(ex, y_zero, 0.0, f_zero, J1);
  odelia::ode::fd_jacobian(ex, y_zero, 0.0, f_zero, J2, 1e-6, 1e-10);
  check(std::abs(J1[0] - 1.0) < 1e-4, "fd_jacobian's default floor gives the unit slope at y = 0");
  check(std::abs(J2[0] - 1.0) > 0.5, "where a floor of 1e-10 loses it (test is not vacuous)");

  // A component on the upper edge of its domain is perturbed downwards when
  // the right-hand side refuses the upward point.
  auto ceiling = [](double, const State& y, State& dydt) {
    if (y[0] > 1.0) odelia::util::stop_domain("over 1");
    dydt[0] = -(y[0] - 1.0);
  };
  CallbackSystem ce(ceiling, State{1.0}, 0.0);
  State y_one{1.0}, f_one(1), Jc;
  odelia::ode::derivs(ce, y_one, f_one, 0.0);
  odelia::ode::fd_jacobian(ce, y_one, 0.0, f_one, Jc);
  check(std::abs(Jc[0] + 1.0) < 1e-4, "fd_jacobian falls back to a downward perturbation on the boundary");
  CallbackSystem ce2(ceiling, State{1.0}, 0.0, CallbackSystem::jac_type(), CallbackSystem::valid_type(), true);
  odelia::ode::Solver<CallbackSystem> sce(ce2, tight_control(), Method::rodas);
  bool threw = false;
  try { sce.advance_adaptive(std::vector<double>{0.0, 1.0}); } catch (const std::runtime_error&) { threw = true; }
  check(!threw && std::abs(sce.state()[0] - 1.0) < 1e-8,
        "so RODAS with a finite-difference Jacobian starts from the boundary");
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

  // A compiled system -- the way plant and leaf use the solver -- under every
  // stepper, from plain C++ with no R anywhere: the method is a constructor
  // argument and nothing else about the system changes.
  {
    LorenzSystem<double> rk(SIGMA, RHO, BETA), dp(SIGMA, RHO, BETA), ro(SIGMA, RHO, BETA);
    odelia::ode::Solver<LorenzSystem<double>> s_rk(rk, tight_control(), Method::rkck);
    odelia::ode::Solver<LorenzSystem<double>> s_dp(dp, tight_control(), Method::dopri);
    odelia::ode::Solver<LorenzSystem<double>> s_ro(ro, tight_control(), Method::rodas);
    s_rk.advance_adaptive(times);
    s_dp.advance_adaptive(times);
    s_ro.advance_adaptive(times);
    check(max_abs_diff(s_dp.state(), s_rk.state()) < 1e-5,
          "a compiled system steps under Dormand-Prince from C++");
    check(max_abs_diff(s_ro.state(), s_rk.state()) < 1e-5,
          "and under RODAS, through its rebind_from() and AD Jacobian");
    std::printf("       (compiled Lorenz, t = 2 at 1e-10: rkck %zu, dopri %zu, rodas %zu steps)\n",
                s_rk.times().size() - 1, s_dp.times().size() - 1, s_ro.times().size() - 1);
  }

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


// --- The step-size rule is a switch (#64) -------------------------------------

void test_controller_switch() {
  using odelia::ode::Controller;
  using odelia::ode::Method;
  const State y0{1.0, 1.0, 1.0};
  const std::vector<double> times{0.0, 20.0};
  odelia::ode::OdeControl control(1e-6, 1e-6, 1.0, 0.0, 1e-12, 1e9, 1e-6);

  // The default is the gsl rule, and saying so changes nothing: bit-identical.
  CallbackSystem a(lorenz_rhs, y0, 0.0);
  odelia::ode::Solver<CallbackSystem> sa(a, control, Method::dopri);
  sa.advance_adaptive(times);
  odelia::ode::OdeControl explicit_gsl = control;
  explicit_gsl.set_controller(Controller::gsl);
  CallbackSystem b(lorenz_rhs, y0, 0.0);
  odelia::ode::Solver<CallbackSystem> sb(b, explicit_gsl, Method::dopri);
  sb.advance_adaptive(times);
  check(control.get_controller() == Controller::gsl, "the default rule is gsl");
  check(sa.state() == sb.state() && sa.times() == sb.times(),
        "and naming it changes nothing, bit for bit");

  // Hairer's rule reaches the same answer (over a horizon short enough that
  // chaos does not separate two correct step sequences) with fewer rejected
  // attempts over a long one.
  odelia::ode::OdeControl hairer = control;
  hairer.set_controller(Controller::hairer);
  {
    CallbackSystem g(lorenz_rhs, y0, 0.0), hh(lorenz_rhs, y0, 0.0);
    odelia::ode::Solver<CallbackSystem> sg(g, control, Method::dopri), sh(hh, hairer, Method::dopri);
    sg.advance_adaptive(std::vector<double>{0.0, 2.0});
    sh.advance_adaptive(std::vector<double>{0.0, 2.0});
    check(max_abs_diff(sh.state(), sg.state()) < 1e-3,
          "the hairer rule integrates Lorenz to the same place at 1e-6");
  }
  CallbackSystem c(lorenz_rhs, y0, 0.0);
  odelia::ode::Solver<CallbackSystem> sc(c, hairer, Method::dopri);
  sc.advance_adaptive(times);
  check(sc.get_n_rejections() < sa.get_n_rejections() * 0.8,
        "with at least a fifth fewer rejected attempts");
  check(sc.get_system_ref().n_rhs <= sa.get_system_ref().n_rhs,
        "and no more evaluations");
  std::printf("       (gsl %zu steps + %zu rejections, hairer %zu + %zu)\n",
              sa.times().size() - 1, sa.get_n_rejections(),
              sc.times().size() - 1, sc.get_n_rejections());

  // The two validity paths are shared: a non-finite estimate and a domain
  // refusal are rejections under either rule.
  auto logistic = [](double, const State& y, State& dydt) {
    if (y[0] < 0.0 || y[0] > 1.0) odelia::util::stop_domain("outside [0, 1]");
    dydt[0] = 50.0 * y[0] * (1.0 - y[0]);
  };
  odelia::ode::OdeControl loose = loose_control();
  loose.set_controller(Controller::hairer);
  CallbackSystem d(logistic, State{0.5}, 0.0);
  odelia::ode::Solver<CallbackSystem> sd(d, loose, Method::rkck);
  sd.advance_adaptive(std::vector<double>{0.0, 1.0});
  check(sd.state()[0] >= 0.0 && sd.state()[0] <= 1.0 && sd.get_n_rejections() > 0,
        "a domain refusal is a rejection under the hairer rule too");
}

} // namespace

int main() {
  std::printf("odelia solver core, standalone (no R, no Rcpp)\n");
  test_stop_throws();
  test_interpolator();
  test_driver_slopes();
  test_solver_runs();
  test_control_rejects_nonfinite_error();
  test_solver_refuses_nonfinite_state();
  test_predicate_rejects_out_of_domain_step();
  test_domain_error_becomes_a_rejection();
  test_non_domain_throw_is_not_absorbed();
  test_unreachable_domain_fails_with_a_reason();
  test_supplied_rows_cost_one_statement();
  test_supplied_rows_carry_a_direction();
  test_a_preaccumulated_residual_keeps_its_rows();
  test_two_preaccumulated_solves_do_not_add_up();

  test_pinned_step_domain_error_is_a_rejection();
  test_pinned_step_predicate_is_enforced();
  test_pinned_step_non_domain_throw_is_not_absorbed();
  test_pinned_step_unreachable_domain_fails_with_a_reason();
  test_pinned_step_unchanged_when_nothing_objects();
  test_sweep_of_a_pinned_run_matches_a_difference();
  test_subdivided_pinned_row_is_refused();
  test_sweep_refuses_one_batch_for_both();
  test_implicit_value_leaves_the_callers_adjoint_alone();
  test_sweep_across_an_insertion();
  test_jacobian_is_closed_at_an_adjoint_scalar();
  test_jacobian_hook_is_detected();
  test_callback_rodas_matches_compiled();
  test_callback_call_budget();
  test_callback_domain_error_is_a_rejection();
  test_callback_single_step_driving();
  test_singular_w_is_a_rejection();
  test_dense_output();
  test_dopri();
  test_controller_switch();
  if (failures == 0) {
    std::printf("all checks passed\n");
    return 0;
  }
  std::printf("%d failure(s)\n", failures);
  return 1;
}
