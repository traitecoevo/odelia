# Tests for the steady-state solve + implicit-function-theorem sensitivity
# capability (issue #39). Exercised through small standalone systems with
# closed-form equilibria and sensitivities, compiled on demand via sourceCpp
# (the in-package Lorenz system has no fixed-point attractor).

ensure_ss_runner <- function() {
  if (isTRUE(.odelia_test_cache$ss_loaded)) {
    return(invisible(TRUE))
  }
  ensure_ode_interface_loaded()

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
    PKG_CPPFLAGS = odelia_cppflags(),
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

# Closed-form interior equilibrium and sensitivity for
#   y0' = r*y0*(1 - y0/K) - m*y0 ; y1' = y0^2 - c*y1
logistic_equilibrium <- function(theta) {
  r <- theta[[1]]; K <- theta[[2]]; m <- theta[[3]]; c <- theta[[4]]
  y0 <- K * (1 - m / r)
  c(y0, y0^2 / c)
}
logistic_sensitivity <- function(theta) {
  r <- theta[[1]]; K <- theta[[2]]; m <- theta[[3]]; c <- theta[[4]]
  y0 <- K * (1 - m / r)
  dy0 <- c(K * m / r^2, 1 - m / r, -K / r, 0)
  dy1 <- c(2 * y0 / c * dy0[1:3], -y0^2 / c^2)
  rbind(dy0, dy1, deparse.level = 0)
}

theta <- c(2.0, 1.5, 0.7)
y_guess <- c(0.0, 0.0)

testthat::test_that("Newton converges to the analytic equilibrium", {
  ensure_ss_runner()
  res <- ss_run(theta, y_guess, warmup = FALSE)

  expect_true(res$converged)
  expect_false(res$warmed)
  expect_lt(res$residual_norm, 1e-10)
  expect_gt(res$iterations, 1L) # genuinely nonlinear: not a one-step solve
  expect_equal(res$y, analytic_equilibrium(theta), tolerance = 1e-10)
})

testthat::test_that("IFT sensitivity matches the closed form", {
  ensure_ss_runner()
  res <- ss_run(theta, y_guess, warmup = FALSE)
  expect_equal(res$sensitivity, analytic_sensitivity(theta), tolerance = 1e-8)
  # The rates read a, b, c live, so the self-check sees only FD noise.
  expect_lt(res$param_check, 1e-6)
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

testthat::test_that("Newton reaches a repelling root and says so", {
  ensure_ss_runner()
  theta_rep <- c(2.0, -1.5, 0.7) # b < 0: y0 runs away from a/b
  res <- ss_run(theta_rep, y_guess, warmup = FALSE)

  expect_true(res$converged)
  expect_equal(res$y, analytic_equilibrium(theta_rep), tolerance = 1e-10)
  expect_equal(res$spectral_abscissa, 1.5, tolerance = 1e-9)
  expect_false(res$stable)
  # The implicit function theorem does not care which kind of root it is.
  expect_equal(res$sensitivity, analytic_sensitivity(theta_rep), tolerance = 1e-8)
})

testthat::test_that("sensitivity and stability refuse after a failed solve", {
  ensure_ss_runner()
  expect_error(ss_run(theta, c(50, 50), warmup = FALSE, max_iter = 1L),
               "converged solve")
})

testthat::test_that("RODAS warm-start reaches the same equilibrium from a far guess", {
  ensure_ss_runner()
  far_guess <- c(100.0, -50.0)
  res <- ss_run(theta, far_guess, warmup = TRUE)

  expect_true(res$converged)
  expect_true(res$stable)
  expect_equal(res$y, analytic_equilibrium(theta), tolerance = 1e-9)
  expect_equal(res$sensitivity, analytic_sensitivity(theta), tolerance = 1e-8)
})

theta_log <- c(1.0, 10.0, 0.4, 0.7)
small_guess <- c(1e-3, 0.0)

testthat::test_that("Newton from a small guess lands on the trivial repelling root", {
  ensure_ss_runner()
  res <- ss_run_logistic(theta_log, small_guess, warmup = FALSE)

  expect_true(res$converged)
  expect_false(res$warmed)
  expect_equal(res$y, c(0, 0), tolerance = 1e-10)
  expect_false(res$stable)
  expect_equal(res$spectral_abscissa, theta_log[1] - theta_log[3], tolerance = 1e-9)
})

testthat::test_that("the warm start gets past the repelling root to the attractor", {
  ensure_ss_runner()
  res <- ss_run_logistic(theta_log, small_guess, warmup = TRUE)

  expect_true(res$converged)
  expect_true(res$warmed)
  expect_true(res$stable)
  expect_equal(res$y, logistic_equilibrium(theta_log), tolerance = 1e-9)
  expect_equal(res$spectral_abscissa,
               max(theta_log[3] - theta_log[1], -theta_log[4]), tolerance = 1e-9)
  expect_equal(res$sensitivity, logistic_sensitivity(theta_log), tolerance = 1e-8)
  expect_lt(res$param_check, 1e-6)
})

testthat::test_that("a parameter read through a cached quantity is caught by the self-check", {
  ensure_ss_runner()
  res <- ss_run_cached(c(2.0, 3.0), 0.0)

  expect_true(res$converged)
  expect_equal(res$y, 6.0, tolerance = 1e-10)
  # The contract breach: the true dy*/d(a, b) is (3, 2); AD seeded in place
  # sees none of k = a*b and reports zeros.
  expect_equal(as.vector(res$sensitivity), c(0, 0))
  expect_gt(res$param_check, 0.5)
})

testthat::test_that("the eigenvalue routine handles what a 2x2 never reaches", {
  ensure_ss_runner()
  check <- function(A, re, im, tol = 1e-12) {
    got <- ss_eigenvalues(A)
    o_got <- order(got$re, got$im)
    o_exp <- order(re, im)
    expect_equal(got$re[o_got], re[o_exp], tolerance = tol)
    expect_equal(got$im[o_got], im[o_exp], tolerance = tol)
  }
  # 3x3 with a complex pair: 1, 2 +/- 3i.
  check(matrix(c(2, -3, 0, 3, 2, 0, 0, 0, 1), 3, byrow = TRUE),
        c(2, 2, 1), c(3, -3, 0))
  # Companion matrix of (x-1)(x-2)(x-3)(x-4).
  check(matrix(c(10, -35, 50, -24, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0), 4, byrow = TRUE),
        c(1, 2, 3, 4), c(0, 0, 0, 0), tol = 1e-10)
  # 5x5 block upper triangular: -1 +/- 2i, -3, -0.5 +/- 2i.
  check(matrix(c(-1, 2, 3, 4, 5,
                 -2, -1, 1, 2, 3,
                 0, 0, -3, 7, 1,
                 0, 0, 0, -0.5, 2,
                 0, 0, 0, -2, -0.5), 5, byrow = TRUE),
        c(-1, -1, -3, -0.5, -0.5), c(2, -2, 0, 2, -2))
  # Defective: a Jordan block, eigenvalue 1 three times.
  check(matrix(c(1, 1, 0, 0, 1, 1, 0, 0, 1), 3, byrow = TRUE),
        c(1, 1, 1), c(0, 0, 0), tol = 1e-5) # a triple root is only O(eps^(1/3))
  # Cyclic permutations: the roots of unity, every one on the unit circle, which
  # a QR iteration without exceptional shifts cannot separate.
  check(rbind(cbind(0, diag(3)), c(1, 0, 0, 0)),
        c(1, -1, 0, 0), c(0, 0, 1, -1))
  s3 <- sqrt(3) / 2
  check(rbind(cbind(0, diag(5)), c(1, 0, 0, 0, 0, 0)),
        c(1, -1, 0.5, 0.5, -0.5, -0.5), c(0, 0, s3, -s3, s3, -s3))
})
