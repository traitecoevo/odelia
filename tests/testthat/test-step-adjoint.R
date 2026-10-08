# Tests for Step::step_adjoint against a finite difference of the same step. The
# adjoint records a tape, so the snippet must link against the odelia shared library
# for the XAD Tape symbols (as the leaf/supplied-derivative interfaces do).

compile_step_adjoint_interface <- function() {
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
    #include <type_traits>
    #include <vector>
    #include <odelia/ode_step.hpp>
    #include <examples/lorenz_system.hpp>

    // Every state handed to a System below, in call order, so the stage states
    // step_adjoint rebuilds can be compared against the ones step() stepped through.
    static std::vector<std::vector<double>> stage_log;

    template <typename T>
    struct RecordingLorenz : LorenzSystem<T> {
      using value_type = T;
      RecordingLorenz(T sigma, T R, T b) : LorenzSystem<T>(sigma, R, b) {}

      template <class S2> using rebind = RecordingLorenz<S2>;

      template <class S2>
      rebind<S2> rebind_from() const {
        std::vector<double> p = this->pars();
        rebind<S2> out{S2(p[0]), S2(p[1]), S2(p[2])};
        std::vector<T> ic(3);
        this->ode_initial_state(ic.begin());
        std::vector<double> ic0 = {xad::value(ic[0]), xad::value(ic[1]),
                                   xad::value(ic[2])};
        out.set_initial_state(ic0.begin(), this->ode_t0());
        return out;
      }

      template <typename Iterator>
      Iterator set_ode_state(Iterator it, double time_ = 0.0) {
        stage_log.push_back({xad::value(it[0]), xad::value(it[1]), xad::value(it[2])});
        return LorenzSystem<T>::set_ode_state(it, time_);
      }
    };

    // The six stage states step() steps through, then the six the recording
    // builds, both ascending in stage. One stage_log serves every scalar the
    // template is instantiated at, so the states the recording builds -- at the
    // adjoint scalar -- land in it beside the ones the forward pass built.
    // [[Rcpp::export]]
    Rcpp::List lorenz_stage_states(std::vector<double> pars, double time,
                                   double step_size, std::vector<double> y,
                                   std::vector<double> lambda_out) {
      RecordingLorenz<double> system(pars[0], pars[1], pars[2]);
      odelia::ode::Step<RecordingLorenz<double>> stepper;
      stepper.resize(y.size());

      std::vector<double> dydt_in(y.size()), dydt_out(y.size()), yerr(y.size());
      std::vector<double> y_end(y);

      // y, then the five later stage states, then y_end.
      stage_log.clear();
      odelia::ode::derivs(system, y, dydt_in, time);
      decltype(stepper)::solved_row solved{};
      stepper.step(system, solved, time, step_size, y_end, yerr, dydt_in,
                   dydt_out);
      Rcpp::List forward(6);
      for (int j = 0; j < 6; ++j) forward[j] = stage_log[j];

      // The six stage states of the recording, ascending, then the start state
      // the step puts the system back to.
      stage_log.clear();
      const odelia::ode::adjoint_rows seeds =
        odelia::ode::adjoint_rows::one_row(lambda_out);
      odelia::ode::adjoint_rows swept;
      odelia::ode::adjoint_rows rows(1, system.ad_parameters().size());
      odelia::ode::adjoint_tape<double> step_tape(false);
      odelia::ode::active_system<std::decay_t<decltype(system)>> active{
          system, step_tape};
      stepper.step_adjoint(active, std::as_const(solved), time, step_size, y, seeds, swept, rows);
      const std::vector<double> lambda_in = swept.to_rows()[0];
      const std::vector<double> parameter_adjoint = rows.to_rows()[0];
      Rcpp::List recorded(6);
      for (int j = 0; j < 6; ++j) recorded[j] = stage_log[j];

      return Rcpp::List::create(Rcpp::_["forward"] = forward,
                                Rcpp::_["recorded"] = recorded,
                                Rcpp::_["last"] = stage_log.back(),
                                Rcpp::_["n_logged"] = (int) stage_log.size());
    }

    // One RKCK step of the Lorenz system from y, plus the adjoint of that step.
    // [[Rcpp::export]]
    Rcpp::List lorenz_step_and_adjoint(std::vector<double> pars, double time,
                                       double step_size, std::vector<double> y,
                                       std::vector<double> lambda_out) {
      LorenzSystem<double> system(pars[0], pars[1], pars[2]);
      odelia::ode::Step<LorenzSystem<double>> stepper;
      stepper.resize(y.size());

      std::vector<double> dydt_in(y.size()), dydt_out(y.size()), yerr(y.size());
      std::vector<double> y_end(y);
      odelia::ode::derivs(system, y, dydt_in, time);
      decltype(stepper)::solved_row solved{};
      stepper.step(system, solved, time, step_size, y_end, yerr, dydt_in,
                   dydt_out);

      const odelia::ode::adjoint_rows seeds =
        odelia::ode::adjoint_rows::one_row(lambda_out);
      odelia::ode::adjoint_rows swept;
      odelia::ode::adjoint_rows rows(1, system.ad_parameters().size());
      odelia::ode::adjoint_tape<double> step_tape(false);
      odelia::ode::active_system<std::decay_t<decltype(system)>> active{
          system, step_tape};
      stepper.step_adjoint(active, std::as_const(solved), time, step_size, y, seeds, swept, rows);
      const std::vector<double> lambda_in = swept.to_rows()[0];
      const std::vector<double> parameter_adjoint = rows.to_rows()[0];

      return Rcpp::List::create(Rcpp::_["y_end"] = y_end,
                                Rcpp::_["lambda_in"] = lambda_in);
    }', verbose = FALSE)
}

testthat::test_that("step_adjoint reproduces a finite difference of one Lorenz step", {
  compile_step_adjoint_interface()

  pars <- c(10.0, 28.0, 8.0 / 3.0)
  time <- 0.0
  step_size <- 0.01
  y <- c(1.5, -0.7, 20.0)
  # Non-trivial and non-symmetric, so a dropped transpose would show up.
  lambda_out <- c(0.3, -1.7, 2.1)
  eps <- 1e-6

  step_end <- function(y) {
    lorenz_step_and_adjoint(pars, time, step_size, y, lambda_out)$y_end
  }

  # J[i, j] = d y_end[i] / d y[j], by central difference of the same step.
  jacobian <- sapply(seq_along(y), function(j) {
    up <- y; up[j] <- up[j] + eps
    down <- y; down[j] <- down[j] - eps
    (step_end(up) - step_end(down)) / (2 * eps)
  })

  expected <- as.vector(t(jacobian) %*% lambda_out)
  got <- lorenz_step_and_adjoint(pars, time, step_size, y, lambda_out)$lambda_in

  expect_equal(length(got), length(y))
  expect_equal(got, expected, tolerance = 1e-7)
  # The transpose is load-bearing: J is not symmetric here.
  expect_false(isTRUE(all.equal(expected, as.vector(jacobian %*% lambda_out),
                                tolerance = 1e-3)))
})

testthat::test_that("the recording builds each stage state bit-identically to step()", {
  compile_step_adjoint_interface()

  r <- lorenz_stage_states(c(10.0, 28.0, 8.0 / 3.0), 0.0, 0.01,
                           c(1.5, -0.7, 20.0), c(0.3, -1.7, 2.1))

  # One state loaded per stage of the tableau, and nothing after them: the stage
  # states are the recording's own intermediates, and the double system is not
  # walked once the recording is taken.
  expect_identical(r$n_logged, 6L)
  # What the sweep transposes has to be the step the forward pass took, term for
  # term: b21 * h * k1 and h * (b21 * k1) round differently, and one function
  # builds both sides so they cannot come apart.
  for (j in seq_len(6)) {
    expect_identical(r$recorded[[j]], r$forward[[j]],
                     info = paste("stage", j))
  }

  # The last state loaded is the sixth stage's: nothing follows the stages.
  expect_identical(r$last, r$recorded[[6]])
})

testthat::test_that("the adjoint scalar is named at namespace scope", {
  withr::local_envvar(PKG_CPPFLAGS = odelia_cppflags())

  # Reached through the interface header alone, with no Solver named and none
  # instantiated: the assertions are on types only.
  Rcpp::sourceCpp(code = '
    // [[Rcpp::plugins(cpp20)]]
    #include <Rcpp.h>
    #include <type_traits>
    #include <odelia/ode_interface.hpp>

    static_assert(std::is_same_v<odelia::ode::active_scalar<>,
                                 xad::adj<double>::active_type>);
    static_assert(std::is_same_v<odelia::ode::active_scalar<float>,
                                 xad::adj<float>::active_type>);

    // [[Rcpp::export]]
    bool active_scalar_named_without_solver() { return true; }
  ')
  expect_true(active_scalar_named_without_solver())
})

# Compile step_adjoint on a minimal System, with `rebind_body` as its rebind_from()
# declaration. Empty gives a System with no hook; one returning Decay gives a hook
# that hands back the wrong scalar. Both must be refused at the call.
compile_step_adjoint_on <- function(rebind_body) {
  withr::local_envvar(PKG_CPPFLAGS = odelia_cppflags())
  Rcpp::sourceCpp(code = sprintf('
    // [[Rcpp::plugins(cpp20)]]
    #include <Rcpp.h>
    #include <type_traits>
    #include <vector>
    #include <odelia/ode_step.hpp>

    struct Decay {
      using value_type = double;
      size_t ode_size() const { return 1; }
      template <typename It> It set_ode_state(It it, double = 0.0) { y = it[0]; return it + 1; }
      template <typename It> It ode_state(It it) const { it[0] = y; return it + 1; }
      template <typename It> It ode_rates(It it) const { it[0] = -y; return it + 1; }
      %s
      double y = 1.0;
    };

    // [[Rcpp::export]]
    std::vector<double> decay_step_adjoint() {
      Decay system;
      odelia::ode::Step<Decay> stepper;
      stepper.resize(1);
      std::vector<double> y(1, 1.0), lambda_out(1, 1.0);
      const odelia::ode::adjoint_rows seeds =
        odelia::ode::adjoint_rows::one_row(lambda_out);
      odelia::ode::adjoint_rows swept;
      odelia::ode::adjoint_rows rows(1, system.ad_parameters().size());
      odelia::ode::adjoint_tape<double> step_tape(false);
      odelia::ode::active_system<std::decay_t<decltype(system)>> active{
          system, step_tape};
      const decltype(stepper)::solved_row solved{};
      stepper.step_adjoint(active, solved, 0.0, 0.1, y, seeds, swept, rows);
      return swept.to_rows()[0];
    }
  ', rebind_body))
}

testthat::test_that("step_adjoint refuses a System that cannot rebind to the adjoint scalar", {
  # No rebind_from() at all.
  expect_error(compile_step_adjoint_on(""))

  # Present, but hands back a System on the wrong scalar.
  expect_error(compile_step_adjoint_on(
    "template <class U> Decay rebind_from() const { return *this; }"))

  # A System that does rebind is stepped, and its adjoint is checked against a
  # finite difference, by the Lorenz tests above.
})
