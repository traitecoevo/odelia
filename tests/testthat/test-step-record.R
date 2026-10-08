# Tests that the solver records the step size each accepted step took, and that a
# replay driven by those sizes reproduces the adaptive run bitwise where a replay
# driven by the recorded times does not.

compile_step_record_interface <- function() {
  ensure_ode_interface_loaded()

  include_dir <- odelia_include_dir()
  odelia_so <- .odelia_test_cache$odelia_so
  pkg_libs <- if (is.character(odelia_so) &&
                  length(odelia_so) == 1 &&
                  !is.na(odelia_so) &&
                  nzchar(odelia_so) &&
                  file.exists(odelia_so)) {
    shQuote(normalizePath(odelia_so, winslash = "/", mustWork = FALSE))
  } else {
    Sys.getenv("PKG_LIBS", unset = "")
  }
  withr::local_envvar(
    PKG_CPPFLAGS = odelia_cppflags(include_dir),
    PKG_LIBS = pkg_libs
  )
  Rcpp::sourceCpp(code = '
    // [[Rcpp::plugins(cpp20)]]
    #include <Rcpp.h>
    #include <vector>
    #include <odelia/ode_solver.hpp>
    #include <examples/lorenz_system.hpp>

    using namespace odelia;
    using LorenzD = LorenzSystem<double>;

    static ode::Solver<LorenzD> make_solver(Rcpp::NumericVector y0, double t0) {
      LorenzD sys(10.0, 28.0, 8.0 / 3.0);
      std::vector<double> y(y0.begin(), y0.end());
      sys.set_initial_state(y.begin(), t0);
      sys.reset();
      ode::Solver<LorenzD> solver(sys, ode::OdeControl());
      solver.set_collect(false);
      return solver;
    }

    // Adaptive run over [t0, t1], then the same trajectory replayed twice: once
    // over the recorded step sizes, once over the recorded times.
    // [[Rcpp::export]]
    Rcpp::List step_record_replays(Rcpp::NumericVector y0, double t0, double t1) {
      ode::Solver<LorenzD> solver = make_solver(y0, t0);
      solver.advance_adaptive({t0, t1});
      const std::vector<double> times = solver.times();
      const std::vector<double> step_sizes = solver.step_sizes();
      const std::vector<double> y_adaptive = solver.state();

      std::vector<ode::instruction> recording;
      recording.reserve(times.size());
      for (size_t i = 0; i < times.size(); ++i) {
        recording.push_back({times[i], step_sizes[i]});
      }
      ode::Solver<LorenzD> by_steps = make_solver(y0, t0);
      by_steps.advance_recorded(recording);

      ode::Solver<LorenzD> by_times = make_solver(y0, t0);
      by_times.advance_fixed(times);

      return Rcpp::List::create(
        Rcpp::Named("times") = Rcpp::wrap(times),
        Rcpp::Named("step_sizes") = Rcpp::wrap(step_sizes),
        Rcpp::Named("y_adaptive") = Rcpp::wrap(y_adaptive),
        Rcpp::Named("y_by_steps") = Rcpp::wrap(by_steps.state()),
        Rcpp::Named("y_by_times") = Rcpp::wrap(by_times.state()),
        Rcpp::Named("time_by_steps") = by_steps.time(),
        Rcpp::Named("time_adaptive") = solver.time());
    }
  ', verbose = FALSE)
}

bits <- function(x) vapply(x, function(v) paste(writeBin(v, raw()), collapse = ""),
                           character(1))

test_that("the solver records a step size beside each recorded time", {
  compile_step_record_interface()
  # t of order 100 against h of order 0.01, where fl(fl(t + h) - t) != h bites.
  res <- step_record_replays(c(1, 1, 1), 100, 101)

  expect_equal(length(res$step_sizes), length(res$times))
  expect_true(is.nan(res$step_sizes[1]))
  expect_true(all(is.finite(res$step_sizes[-1])))
  expect_true(all(res$step_sizes[-1] > 0))

  # sum(h) reaches the final time only approximately; that is why h is recorded.
  reached <- res$times[1] + sum(res$step_sizes[-1])
  expect_equal(reached, 101, tolerance = 1e-12)
  cat(sprintf("steps=%d h range=[%.6g, %.6g] t0+sum(h)-t1=%.3g\n",
              length(res$step_sizes) - 1L, min(res$step_sizes[-1]),
              max(res$step_sizes[-1]), reached - 101))
})

test_that("a replay over the recorded step sizes reproduces the adaptive run bitwise", {
  compile_step_record_interface()
  res <- step_record_replays(c(1, 1, 1), 100, 101)

  expect_identical(bits(res$y_by_steps), bits(res$y_adaptive))

  # And the time-driven replay does not: differencing the recorded times cannot
  # recover the step sizes the adaptive run took.
  differing <- sum(bits(res$y_by_times) != bits(res$y_adaptive))
  worst <- max(abs(res$y_by_times - res$y_adaptive))
  cat(sprintf("time-driven replay: %d of %d components differ, worst |diff|=%.6g\n",
              differing, length(res$y_adaptive), worst))
  expect_gt(differing, 0)
})

