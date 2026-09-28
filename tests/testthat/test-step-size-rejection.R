# Tests for the adaptive controller's step decision on a non-finite error estimate.
# The per-component error ratios are reduced with std::max, which does not propagate
# NaN: a non-finite ratio followed by a finite one is dropped, so the largest measured
# error comes from whichever components stayed finite and the step is accepted -- and
# grown, if those are small. With no finite ratio after it rmax is NaN, both the
# too-large and too-small tests are false, and the step is accepted unchanged. Either
# way a non-finite state enters the trajectory.
#
# Driven from a System through the solver's own step loop, so what is measured is the
# decision the solver acts on. The solver template needs the odelia shared library for
# the XAD Tape symbols, so this snippet links against it rather than compiling
# header-only.

compile_step_size_rejection_interface <- function() {
  ensure_ode_interface_loaded()

  include_dir <- odelia_include_dir()
  odelia_so <- .odelia_test_cache$odelia_so
  withr::local_envvar(
    PKG_CPPFLAGS = odelia_cppflags(include_dir),
    PKG_LIBS = shQuote(normalizePath(odelia_so, winslash = "/", mustWork = TRUE))
  )
  Rcpp::sourceCpp(code = '
    // [[Rcpp::plugins(cpp20)]]
    #include <Rcpp.h>
    #include <vector>
    #include <cmath>
    #include <limits>
    #include <odelia/ode_solver.hpp>

    // A one-state system whose rate goes non-finite once the state passes a
    // threshold, so the adaptive controller meets a non-finite error estimate.
    struct DivergingSystem {
      using value_type = double;
      double y = 1.0, dydt = 1.0, time = 0.0, t0 = 0.0, threshold = 1.5;

      size_t ode_size() const { return 1; }
      double ode_time() const { return time; }
      double ode_t0() const { return t0; }

      template <typename Iterator>
      Iterator set_ode_state(Iterator it, double time_) {
        y = *it++;
        time = time_;
        dydt = (y > threshold) ? std::numeric_limits<double>::quiet_NaN() : 1.0;
        return it;
      }
      template <typename Iterator>
      Iterator set_initial_state(Iterator it, double t0_ = 0.0) {
        t0 = t0_;
        y = *it++;
        return it;
      }
      template <typename Iterator>
      Iterator ode_state(Iterator it) const { *it++ = y; return it; }
      template <typename Iterator>
      Iterator ode_rates(Iterator it) const { *it++ = dydt; return it; }
      void reset() { y = 1.0; time = t0; dydt = 1.0; }
    };

    // Advance the diverging system past its threshold and report the state reached.
    // [[Rcpp::export]]
    std::vector<double> advance_past_non_finite(double t_end) {
      odelia::ode::OdeControl control;
      odelia::ode::Solver<DivergingSystem> solver(DivergingSystem(), control);
      solver.set_state(std::vector<double>{1.0}, 0.0);
      solver.advance_adaptive(std::vector<double>{0.0, t_end});
      return solver.state();
    }', verbose = FALSE)
}

testthat::test_that("a non-finite step-size decision is rejected", {
  compile_step_size_rejection_interface()

  # Every attempt past the threshold is rejected and shrunk by the largest permitted
  # factor, so the step loop runs down to the accuracy limit and stops. Accepting
  # instead would commit a non-finite state to the trajectory, which a replay pinned
  # to that grid then differentiates.
  testthat::expect_error(advance_past_non_finite(2.0),
                         "Cannot achieve the desired accuracy")

  # Below the threshold the rates stay finite and the same run advances normally, so
  # the rejection costs the finite path nothing.
  testthat::expect_equal(advance_past_non_finite(0.25), 1.25, tolerance = 1e-10)
})
