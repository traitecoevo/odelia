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

} // namespace

int main() {
  std::printf("odelia solver core, standalone (no R, no Rcpp)\n");
  test_stop_throws();
  test_interpolator();
  test_solver_runs();
  if (failures == 0) {
    std::printf("all checks passed\n");
    return 0;
  }
  std::printf("%d failure(s)\n", failures);
  return 1;
}
