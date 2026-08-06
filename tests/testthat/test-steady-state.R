# Tests for the steady-state solve + implicit-function-theorem sensitivity
# capability (issue #39). Exercised through a small standalone system with a
# closed-form equilibrium and sensitivity, compiled on demand via sourceCpp
# (the in-package Lorenz system has no fixed-point attractor).

ensure_ss_runner <- function() {
  if (isTRUE(.odelia_test_cache$ss_loaded)) {
    return(invisible(TRUE))
  }
  ensure_ode_interface_loaded()

  include_dir <- dirname(dirname(resolve_test_path(
    "include/odelia/ode_solver.hpp", "inst/include/odelia/ode_solver.hpp")))
  runner_cpp <- resolve_test_path(
    "tests/testthat/steady_state_runner.cpp",
    "tests/testthat/steady_state_runner.cpp")

  odelia_so <- .odelia_test_cache$odelia_so
  pkg_libs <- if (is.character(odelia_so) && length(odelia_so) == 1 &&
                  !is.na(odelia_so) && nzchar(odelia_so) &&
                  file.exists(odelia_so)) {
    shQuote(normalizePath(odelia_so, winslash = "/", mustWork = FALSE))
  } else {
    Sys.getenv("PKG_LIBS", unset = "")
  }
  withr::local_envvar(
    PKG_CPPFLAGS = paste0("-I", include_dir),
    PKG_LIBS = pkg_libs
  )

  res <- tryCatch({
    Rcpp::sourceCpp(runner_cpp, verbose = FALSE)
    NULL
  }, error = function(e) e)
  if (inherits(res, "error")) {
    if (grepl("active_tape_", conditionMessage(res), fixed = TRUE)) {
      testthat::skip("Steady-state runner symbols unavailable in this load_all session.")
    }
    stop(res)
  }

  .odelia_test_cache$ss_loaded <- TRUE
  invisible(TRUE)
}

# Closed-form equilibrium and sensitivity for
#   y0' = a - b*y0 ; y1' = y0^2 - c*y1
analytic_equilibrium <- function(theta) {
  a <- theta[[1]]; b <- theta[[2]]; c <- theta[[3]]
  c(a / b, a^2 / (b^2 * c))
}
analytic_sensitivity <- function(theta) {
  a <- theta[[1]]; b <- theta[[2]]; c <- theta[[3]]
  matrix(c(
    1 / b,             -a / b^2,             0,
    2 * a / (b^2 * c), -2 * a^2 / (b^3 * c), -a^2 / (b^2 * c^2)
  ), nrow = 2, byrow = TRUE)
}

theta <- c(2.0, 1.5, 0.7)
y_guess <- c(0.0, 0.0)

testthat::test_that("Newton converges to the analytic equilibrium", {
  ensure_ss_runner()
  res <- ss_run(theta, y_guess, warmup = FALSE)

  expect_true(res$converged)
  expect_lt(res$residual_norm, 1e-10)
  expect_gt(res$iterations, 1L) # genuinely nonlinear: not a one-step solve
  expect_equal(res$y, analytic_equilibrium(theta), tolerance = 1e-10)
})

testthat::test_that("IFT sensitivity matches the closed form", {
  ensure_ss_runner()
  res <- ss_run(theta, y_guess, warmup = FALSE)
  expect_equal(res$sensitivity, analytic_sensitivity(theta), tolerance = 1e-8)
})

testthat::test_that("IFT sensitivity matches finite differences", {
  ensure_ss_runner()
  res <- ss_run(theta, y_guess, warmup = FALSE)

  # Central finite differences of the (independently re-solved) equilibrium
  # w.r.t. each parameter -- the validation the issue asks for.
  fd <- matrix(0, nrow = 2, ncol = 3)
  for (j in seq_len(3)) {
    h <- 1e-6 * max(abs(theta[j]), 1)
    tp <- theta; tp[j] <- tp[j] + h
    tm <- theta; tm[j] <- tm[j] - h
    fd[, j] <- (ss_equilibrium(tp, y_guess) - ss_equilibrium(tm, y_guess)) / (2 * h)
  }
  expect_equal(res$sensitivity, fd, tolerance = 1e-6)
})

testthat::test_that("stability is read off the factored Jacobian", {
  ensure_ss_runner()
  res <- ss_run(theta, y_guess, warmup = FALSE)

  # Eigenvalues of df/dy are exactly -b and -c.
  expect_equal(sort(res$eig_re), sort(c(-theta[2], -theta[3])), tolerance = 1e-9)
  expect_true(all(abs(res$eig_im) < 1e-9))
  expect_equal(res$spectral_abscissa, -min(theta[2], theta[3]), tolerance = 1e-9)
  expect_true(res$stable)
  # Autonomous system: no explicit time dependence.
  expect_lt(res$time_dependence, 1e-6)
})

testthat::test_that("RODAS warm-start reaches the same equilibrium from a far guess", {
  ensure_ss_runner()
  far_guess <- c(100.0, -50.0)
  res <- ss_run(theta, far_guess, warmup = TRUE)

  expect_true(res$converged)
  expect_equal(res$y, analytic_equilibrium(theta), tolerance = 1e-9)
  expect_equal(res$sensitivity, analytic_sensitivity(theta), tolerance = 1e-8)
})
