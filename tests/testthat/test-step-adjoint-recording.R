# Tests for the recording one step is transposed by -- all six Runge-Kutta stages
# and the combination closing them, taken on a System rebound to the adjoint
# scalar for that recording -- and for Solver::solve_adjoint over the recorded
# steps. The adjoint records a tape, so the snippet must link against the odelia
# shared library for the XAD Tape symbols.

# Predator and prey, with the four coefficients declared as the parameters a
# recording carries. `dndt` and `dpdt` are DERIVED: written from an expression
# rather than handed a fresh value, which is what a System carried between
# recordings gets wrong.
lv_system <- '
  static int lv_assigns = 0;

  template <typename T>
  class LotkaVolterra {
  public:
    using value_type = T;
    LotkaVolterra(T a_ = T(0), T b_ = T(0), T c_ = T(0), T d_ = T(0))
      : a(a_), b(b_), c(c_), d(d_) { reset(); }

    template <typename> friend class LotkaVolterra;

    // The one map, reached through rebind_from below.
    template <class S2>
    void assign_from(const LotkaVolterra<S2>& src) {
      ++lv_assigns;
      a = T(odelia::util::to_passive(src.a));
      b = T(odelia::util::to_passive(src.b));
      c = T(odelia::util::to_passive(src.c));
      d = T(odelia::util::to_passive(src.d));
      n = T(odelia::util::to_passive(src.n));
      p = T(odelia::util::to_passive(src.p));
      time = src.time;
      compute_rates();
    }

    template <class S2>
    LotkaVolterra<S2> rebind_from() const {
      LotkaVolterra<S2> out;
      out.assign_from(*this);
      return out;
    }

    std::vector<T*> ad_parameters() { return {&a, &b, &c, &d}; }

    // Everything carrying the scalar. The rates are among them: a recording
    // writes them, and writing leaves whatever slot the value already held.
    template <class F>
    void for_each_active(F&& f) {
      f(a); f(b); f(c); f(d);
      f(n); f(p);
      f(dndt); f(dpdt);
    }

    size_t ode_size() const { return 2; }
    double ode_time() const { return time; }
    void reset() { n = 10.0; p = 5.0; time = 0.0; compute_rates(); }

    template <typename It> It set_ode_state(It it, double time_) {
      time = time_;
      n = *it++;
      p = *it++;
      compute_rates();
      return it;
    }
    // This width never moves, so being the shape a recorded time implies is
    // nothing, and the load is the ordinary one.
    void set_recorded_state(const std::vector<T>& y, double time_) {
      set_ode_state(y.begin(), time_);
    }
    void compute_rates() {
      ++rate_calls;
      const T flux = b * n * p;
      dndt = a * n - flux;
      dpdt = c * flux - d * p;
    }
    template <typename It> It ode_state(It it) const {
      *it++ = n; *it++ = p; return it;
    }
    template <typename It> It ode_rates(It it) const {
      *it++ = dndt; *it++ = dpdt; return it;
    }

    int rate_calls = 0;

  private:
    T a, b, c, d;
    T n = 0, p = 0, dndt = 0, dpdt = 0;
    double time = 0;
  };
'

compile_recording_interface <- function() {
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
  Rcpp::sourceCpp(code = paste0('
    // [[Rcpp::plugins(cpp20)]]
    #include <Rcpp.h>
    #include <vector>
    #include <odelia/ode_solver.hpp>
    #include <examples/lorenz_system.hpp>
    ', lv_system, '

    // Both systems in this file are rebound to the adjoint scalar by the
    // recording, so both have to declare the hook it rebinds them with.
    static_assert(odelia::ode::Rebindable<LotkaVolterra<double>,
                                          odelia::ode::active_scalar<double> >);
    static_assert(odelia::ode::Rebindable<LorenzSystem<double>,
                                          odelia::ode::active_scalar<double> >);

    // [[Rcpp::export]]
    bool both_systems_rebind() { return true; }

    // One RKCK step and its adjoint, plus what the step cost and how many times
    // a System was copied for it.
    // [[Rcpp::export]]
    Rcpp::List lv_step_and_adjoint(std::vector<double> pars, double time,
                                   double step_size, std::vector<double> y,
                                   std::vector<double> lambda_out) {
      LotkaVolterra<double> system(pars[0], pars[1], pars[2], pars[3]);
      odelia::ode::Step<LotkaVolterra<double> > stepper;
      stepper.resize(y.size());

      std::vector<double> dydt_in(y.size()), dydt_out(y.size()), yerr(y.size());
      std::vector<double> y_end(y);
      odelia::ode::derivs(system, y, dydt_in, time);
      decltype(stepper)::solved_row solved{};
      stepper.step(system, solved, time, step_size, y_end, yerr, dydt_in,
                   dydt_out);

      LotkaVolterra<double> adj(pars[0], pars[1], pars[2], pars[3]);
      adj.rate_calls = 0;
      lv_assigns = 0;
      const odelia::ode::adjoint_rows seeds =
        odelia::ode::adjoint_rows::one_row(lambda_out);
      odelia::ode::adjoint_rows swept;
      odelia::ode::adjoint_rows rows(1, adj.ad_parameters().size());
      odelia::ode::adjoint_tape<double> step_tape(false);
      odelia::ode::active_system<LotkaVolterra<double>> active{adj, step_tape};
      stepper.step_adjoint(active, std::as_const(solved), time, step_size, y, seeds, swept, rows);
      const std::vector<double> lambda_in = swept.to_rows()[0];
      const std::vector<double> parameter_adjoint = rows.to_rows()[0];

      return Rcpp::List::create(Rcpp::_["y_end"] = y_end,
                                Rcpp::_["lambda_in"] = lambda_in,
                                Rcpp::_["parameter_adjoint"] = parameter_adjoint,
                                Rcpp::_["rate_calls"] = adj.rate_calls,
                                Rcpp::_["assigns"] = lv_assigns,
                                Rcpp::_["recorded_rates"] =
                                  (int) stepper.recorded_rates);
    }

    // The recorded run, the states it passed through, and the adjoint swept back
    // over them.
    // [[Rcpp::export]]
    Rcpp::List lv_solve_adjoint(std::vector<double> pars, std::vector<double> y0,
                                double t_end, std::vector<double> lambda_end) {
      odelia::ode::OdeControl ctl;
      LotkaVolterra<double> system(pars[0], pars[1], pars[2], pars[3]);
      odelia::ode::Solver<LotkaVolterra<double> > solver(system, ctl);
      solver.set_state(y0, 0.0);
      solver.advance_adaptive({0.0, t_end});
      const std::vector<odelia::ode::instruction> recording = solver.schedule();

      // Replays the recorded sizes, keeping its states, so the sweep reads the
      // record the replay itself kept -- rather than copying a whole System per
      // step out of the history to pull one state vector back out of it.
      odelia::ode::Solver<LotkaVolterra<double> > replay(system, ctl);
      replay.set_keep_states(true);
      replay.set_state(y0, 0.0);
      replay.advance_recorded(recording);
      const auto rec = replay.recording();

      odelia::ode::adjoint_rows lambda = odelia::ode::adjoint_rows::one_row(lambda_end);
      odelia::ode::adjoint_rows rows(1, system.ad_parameters().size());
      replay.clear_recorded_rates();
      replay.solve_adjoint(lambda, rows);
      const std::vector<double> parameter_adjoint = rows.to_rows()[0];

      return Rcpp::List::create(Rcpp::_["n_steps"] = (int) recording.size(),
                                Rcpp::_["y_end"] = rec.back().state,
                                Rcpp::_["lambda"] = lambda.to_rows()[0],
                                Rcpp::_["parameter_adjoint"] = parameter_adjoint,
                                Rcpp::_["recorded_rates"] =
                                  (int) replay.recorded_rates());
    }

    // The same sweep taken as two segments meeting at `split`, which is what a
    // caller whose System changes width between steps does.
    // [[Rcpp::export]]
    std::vector<double> lv_solve_adjoint_segments(std::vector<double> pars,
                                                  std::vector<double> y0,
                                                  double t_end,
                                                  std::vector<double> lambda_end,
                                                  int split) {
      odelia::ode::OdeControl ctl;
      LotkaVolterra<double> system(pars[0], pars[1], pars[2], pars[3]);
      odelia::ode::Solver<LotkaVolterra<double> > solver(system, ctl);
      solver.set_state(y0, 0.0);
      solver.advance_adaptive({0.0, t_end});
      const std::vector<odelia::ode::instruction> recording = solver.schedule();

      // Keeps its states, so the sweep reads the record this replay kept.
      odelia::ode::Solver<LotkaVolterra<double> > replay(system, ctl);
      replay.set_keep_states(true);
      replay.set_state(y0, 0.0);
      replay.advance_recorded(recording);
      const auto rec = replay.recording();

      odelia::ode::adjoint_rows lambda = odelia::ode::adjoint_rows::one_row(lambda_end);
      odelia::ode::adjoint_rows rows(1, system.ad_parameters().size());
      replay.solve_adjoint(lambda, rows, (size_t) split, rec.size() - 1);
      replay.solve_adjoint(lambda, rows, 0, (size_t) split);
      return lambda.to_rows()[0];
    }

    // The schedule `schedule_pars` resolved, run forward at `run_pars`. The
    // sweep treats the recorded step sizes as constant, so a difference in the
    // coefficients has to hold them at the values the sweep replayed; letting
    // them move would difference the schedule as well and referee something the
    // sweep does not compute.
    // [[Rcpp::export]]
    std::vector<double> lv_replay_at(std::vector<double> schedule_pars,
                                     std::vector<double> run_pars,
                                     std::vector<double> y0, double t_end) {
      odelia::ode::OdeControl ctl;
      LotkaVolterra<double> schedule_system(schedule_pars[0], schedule_pars[1],
                                            schedule_pars[2], schedule_pars[3]);
      odelia::ode::Solver<LotkaVolterra<double> > solver(schedule_system, ctl);
      solver.set_state(y0, 0.0);
      solver.advance_adaptive({0.0, t_end});
      const std::vector<odelia::ode::instruction> recording = solver.schedule();

      LotkaVolterra<double> run_system(run_pars[0], run_pars[1], run_pars[2],
                                       run_pars[3]);
      odelia::ode::Solver<LotkaVolterra<double> > replay(run_system, ctl);
      replay.set_state(y0, 0.0);
      replay.advance_recorded(recording);
      return replay.state();
    }

    // The same recorded schedule run forward from a given start state, so a
    // finite difference of the whole run reads the identical step sizes.
    // [[Rcpp::export]]
    std::vector<double> lv_replay(std::vector<double> pars,
                                  std::vector<double> y0, double t_end,
                                  std::vector<double> y_start) {
      odelia::ode::OdeControl ctl;
      LotkaVolterra<double> system(pars[0], pars[1], pars[2], pars[3]);
      odelia::ode::Solver<LotkaVolterra<double> > solver(system, ctl);
      solver.set_state(y0, 0.0);
      solver.advance_adaptive({0.0, t_end});
      const std::vector<odelia::ode::instruction> recording = solver.schedule();

      odelia::ode::Solver<LotkaVolterra<double> > replay(system, ctl);
      replay.set_state(y_start, 0.0);
      replay.advance_recorded(recording);
      return replay.state();
    }

    // The same one step, run forward at perturbed coefficients, so a difference
    // of it referees the parameter half of the recording.
    // [[Rcpp::export]]
    std::vector<double> lv_step_end(std::vector<double> pars, double time,
                                    double step_size, std::vector<double> y) {
      LotkaVolterra<double> system(pars[0], pars[1], pars[2], pars[3]);
      odelia::ode::Step<LotkaVolterra<double> > stepper;
      stepper.resize(y.size());
      std::vector<double> dydt_in(y.size()), dydt_out(y.size()), yerr(y.size());
      std::vector<double> y_end(y);
      odelia::ode::derivs(system, y, dydt_in, time);
      decltype(stepper)::solved_row solved{};
      stepper.step(system, solved, time, step_size, y_end, yerr, dydt_in,
                   dydt_out);
      return y_end;
    }
  '), verbose = FALSE)
}

lv_pars <- c(1.1, 0.06, 0.4, 0.9)
lv_y <- c(11.0, 4.0)

testthat::test_that("both systems declare the hook the recording rebinds them with", {
  compile_recording_interface()
  expect_true(both_systems_rebind())
})

testthat::test_that("step_adjoint matches a finite difference of one step", {
  compile_recording_interface()

  step_size <- 0.05
  # Non-symmetric, so a dropped transpose shows up.
  lambda_out <- c(0.7, -1.9)
  eps <- 1e-7

  step_end <- function(y) lv_step_end(lv_pars, 0.0, step_size, y)
  jacobian <- sapply(seq_along(lv_y), function(j) {
    up <- lv_y; up[j] <- up[j] + eps
    down <- lv_y; down[j] <- down[j] - eps
    (step_end(up) - step_end(down)) / (2 * eps)
  })
  expected <- as.vector(t(jacobian) %*% lambda_out)

  got <- lv_step_and_adjoint(lv_pars, 0.0, step_size, lv_y, lambda_out)$lambda_in
  expect_equal(got, expected, tolerance = 1e-6)
  # The transpose is load-bearing: J is not symmetric here.
  expect_false(isTRUE(all.equal(expected, as.vector(jacobian %*% lambda_out),
                                tolerance = 1e-3)))
})

testthat::test_that("the coefficients' rows come back with the state's", {
  compile_recording_interface()

  # The parameters ride in the same recording as the state, so an active System carried
  # into a second recording holds slots that recording has cleared and these
  # rows are what comes back wrong. A stale active System leaves the state rows exact.
  step_size <- 0.05
  lambda_out <- c(0.7, -1.9)
  eps <- 1e-7

  expected <- sapply(seq_along(lv_pars), function(j) {
    up <- lv_pars; up[j] <- up[j] + eps
    down <- lv_pars; down[j] <- down[j] - eps
    sum(((lv_step_end(up, 0.0, step_size, lv_y) -
          lv_step_end(down, 0.0, step_size, lv_y)) / (2 * eps)) * lambda_out)
  })

  got <- lv_step_and_adjoint(lv_pars, 0.0, step_size, lv_y,
                             lambda_out)$parameter_adjoint
  expect_equal(got, expected, tolerance = 1e-6)
  expect_false(any(abs(expected) < 1e-8))
})

testthat::test_that("one component of lambda_out at a time agrees, row by row", {
  compile_recording_interface()

  step_size <- 0.05
  eps <- 1e-7
  for (i in seq_along(lv_y)) {
    seed <- rep(0, length(lv_y))
    seed[i] <- 1
    row <- sapply(seq_along(lv_y), function(j) {
      up <- lv_y; up[j] <- up[j] + eps
      down <- lv_y; down[j] <- down[j] - eps
      (lv_step_end(lv_pars, 0.0, step_size, up)[i] -
       lv_step_end(lv_pars, 0.0, step_size, down)[i]) / (2 * eps)
    })
    got <- lv_step_and_adjoint(lv_pars, 0.0, step_size, lv_y, seed)$lambda_in
    expect_equal(got, row, tolerance = 1e-6,
                 info = paste("seeded component", i))
  }
})

testthat::test_that("a step is one recording of six rates, on one copy of the System", {
  compile_recording_interface()

  r <- lv_step_and_adjoint(lv_pars, 0.0, 0.05, lv_y, c(0.7, -1.9))
  # One recording spans the step, so the System is rebound once for it. Every
  # scalar a recording writes has to arrive holding no tape slot, and a copy
  # taken for the recording is what gives it that.
  expect_identical(r$assigns, 1L)
  # One rate evaluation per stage of the tableau, and the stage states between
  # them are the recording's own intermediates.
  expect_identical(r$recorded_rates, 6L)
  # The double System is not walked. The stages run on the active copy, and a
  # sweep is handed each step's state, so this System is read for its width alone.
  expect_identical(r$rate_calls, 0L)
})

testthat::test_that("solve_adjoint over the recorded steps matches a finite difference of the run", {
  compile_recording_interface()

  t_end <- 1.5
  lambda_end <- c(0.4, -1.3)
  eps <- 1e-6

  r <- lv_solve_adjoint(lv_pars, lv_y, t_end, lambda_end)
  expect_gt(r$n_steps, 3)
  # Six rate evaluations recorded per step swept. The recording's first entry is
  # the state no step reached, so the sweep takes one fewer.
  expect_identical(r$recorded_rates, 6L * (r$n_steps - 1L))

  jacobian <- sapply(seq_along(lv_y), function(j) {
    up <- lv_y; up[j] <- up[j] + eps
    down <- lv_y; down[j] <- down[j] - eps
    (lv_replay(lv_pars, lv_y, t_end, up) -
     lv_replay(lv_pars, lv_y, t_end, down)) / (2 * eps)
  })
  expected <- as.vector(t(jacobian) %*% lambda_end)
  expect_equal(r$lambda, expected, tolerance = 1e-5)
})

testthat::test_that("the coefficients' rows over a whole run match a finite difference of it", {
  compile_recording_interface()

  # The state channel above and the coefficients' channel below fail
  # differently, and only over more than one step. A step composed wrongly
  # leaves the state rows out by a percent or so and the coefficients' rows out
  # by tens of percent, both finite and plausible -- so a suite refereeing the
  # coefficients at one step and the state over a run reports nothing. Both
  # channels are needed over a run.
  t_end <- 1.5
  lambda_end <- c(0.4, -1.3)
  eps <- 1e-6

  got <- lv_solve_adjoint(lv_pars, lv_y, t_end, lambda_end)$parameter_adjoint
  expected <- sapply(seq_along(lv_pars), function(j) {
    up <- lv_pars; up[j] <- up[j] + eps
    down <- lv_pars; down[j] <- down[j] - eps
    sum(((lv_replay_at(lv_pars, up, lv_y, t_end) -
          lv_replay_at(lv_pars, down, lv_y, t_end)) / (2 * eps)) * lambda_end)
  })

  expect_equal(got, expected, tolerance = 1e-5)
  # Every row carries something, so agreement is not two vectors of zeros.
  expect_false(any(abs(expected) < 1e-6))
})

testthat::test_that("the sweep taken as two segments equals the sweep taken whole", {
  compile_recording_interface()

  t_end <- 1.5
  lambda_end <- c(0.4, -1.3)
  whole <- lv_solve_adjoint(lv_pars, lv_y, t_end, lambda_end)

  # Every interior split, so no segment is empty and every boundary is a step
  # boundary. Bit-identity, not a tolerance: a segment sweeps the same steps in
  # the same order and the arithmetic does not change.
  for (split in seq_len(whole$n_steps - 2L)) {
    got <- lv_solve_adjoint_segments(lv_pars, lv_y, t_end, lambda_end, split)
    expect_identical(got, whole$lambda, info = paste("split at", split))
  }
})

testthat::test_that("a segment that is not a range of recorded steps is refused", {
  compile_recording_interface()

  t_end <- 1.5
  lambda_end <- c(0.4, -1.3)
  n <- lv_solve_adjoint(lv_pars, lv_y, t_end, lambda_end)$n_steps
  # A split at either end leaves one of the two segments empty.
  expect_error(lv_solve_adjoint_segments(lv_pars, lv_y, t_end, lambda_end, 0L),
               "not a range of recorded steps")
  expect_error(lv_solve_adjoint_segments(lv_pars, lv_y, t_end, lambda_end,
                                         n - 1L),
               "not a range of recorded steps")
})
