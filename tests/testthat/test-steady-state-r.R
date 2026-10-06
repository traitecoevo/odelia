# The steady-state solve over an R right-hand side (#39): ode_steady_state().
# The solver itself is proven against closed forms through compiled Systems in
# test-steady-state.R; these tests prove the R surface against the same
# systems, with the Jacobian by finite differences or a supplied jacfunc and
# the sensitivity by central differences through func.

# y0' = a - b*y0 ; y1' = y0^2 - c*y1, deSolve shape.
demog_rhs <- function(t, y, p) {
  list(c(p[["a"]] - p[["b"]] * y[[1]], y[[1]]^2 - p[["c"]] * y[[2]]))
}
demog_jac <- function(t, y, p) {
  matrix(c(-p[["b"]], 2 * y[[1]], 0, -p[["c"]]), 2, 2)
}
demog_pars <- c(a = 2.0, b = 1.5, c = 0.7)
demog_equilibrium <- function(p) {
  c(p[["a"]] / p[["b"]], p[["a"]]^2 / (p[["b"]]^2 * p[["c"]]))
}
demog_sensitivity <- function(p) {
  a <- p[["a"]]; b <- p[["b"]]; c <- p[["c"]]
  matrix(c(1 / b, -a / b^2, 0,
           2 * a / (b^2 * c), -2 * a^2 / (b^3 * c), -a^2 / (b^2 * c^2)),
         2, 3, byrow = TRUE)
}

# y0' = r*y0*(1 - y0/K) - m*y0 ; y1' = y0^2 - c*y1: a trivial root at 0.
logistic_rhs <- function(t, y, p) {
  c(p[["r"]] * y[[1]] * (1 - y[[1]] / p[["K"]]) - p[["m"]] * y[[1]],
    y[[1]]^2 - p[["c"]] * y[[2]])
}
logistic_pars <- c(r = 1.0, K = 10.0, m = 0.4, c = 0.7)
logistic_equilibrium <- function(p) {
  y0 <- p[["K"]] * (1 - p[["m"]] / p[["r"]])
  c(y0, y0^2 / p[["c"]])
}
logistic_sensitivity <- function(p) {
  r <- p[["r"]]; K <- p[["K"]]; m <- p[["m"]]; c <- p[["c"]]
  y0 <- K * (1 - m / r)
  dy0 <- c(K * m / r^2, 1 - m / r, -K / r, 0)
  rbind(dy0, c(2 * y0 / c * dy0[1:3], -y0^2 / c^2), deparse.level = 0)
}

testthat::test_that("ode_steady_state reaches the closed-form fixed point, with its stability and sensitivity", {
  eq <- ode_steady_state(demog_rhs, y0 = c(n = 0, m = 0), parms = demog_pars)
  expect_s3_class(eq, "odelia_steady_state")
  expect_true(eq$converged)
  expect_false(eq$warmed)
  expect_gt(eq$iterations, 1)
  expect_lt(eq$residual_norm, 1e-10)
  expect_equal(eq$y, c(n = 2 / 1.5, m = 4 / (1.5^2 * 0.7)), tolerance = 1e-9)
  expect_named(eq$residual, c("n", "m"))
  # The Jacobian is by forward differences through func, so 1e-6-ish.
  expect_equal(dimnames(eq$jacobian), list(c("n", "m"), c("n", "m")))
  expect_equal(unname(eq$jacobian), matrix(c(-1.5, 2 * eq$y[[1]], 0, -0.7), 2, 2),
               tolerance = 1e-5)
  expect_type(eq$eigenvalues, "complex")
  expect_equal(sort(Re(eq$eigenvalues)), c(-1.5, -0.7), tolerance = 1e-5)
  expect_true(all(abs(Im(eq$eigenvalues)) < 1e-9))
  expect_equal(eq$spectral_abscissa, -0.7, tolerance = 1e-5)
  expect_true(eq$stable)
  expect_lt(eq$time_dependence, 1e-6)
  expect_equal(dimnames(eq$sensitivity), list(c("n", "m"), c("a", "b", "c")))
  expect_equal(unname(eq$sensitivity), demog_sensitivity(demog_pars), tolerance = 1e-5)
  expect_named(ode_counts(eq), c("n_rhs", "n_jac", "n_iterations"))
  expect_equal(ode_counts(eq)[["n_iterations"]], eq$iterations)
  expect_output(print(eq), "converged")
})

testthat::test_that("a supplied jacfunc gives the same answer for fewer evaluations", {
  fd <- ode_steady_state(demog_rhs, c(0, 0), demog_pars)
  an <- ode_steady_state(demog_rhs, c(0, 0), demog_pars, jacfunc = demog_jac)
  expect_equal(an$y, fd$y, tolerance = 1e-9)
  expect_equal(an$jacobian, fd$jacobian, tolerance = 1e-5)
  expect_equal(an$sensitivity, fd$sensitivity, tolerance = 1e-5)
  # One Jacobian per iteration plus one at the solution, formed by one call to
  # jacfunc each; by differences each costs n = 2 evaluations of func instead
  # (and an inexact Jacobian may cost Newton an iteration).
  expect_equal(ode_counts(an)[["n_jac"]], an$iterations + 1)
  expect_gte(ode_counts(fd)[["n_rhs"]] - ode_counts(an)[["n_rhs"]], 2 * (an$iterations + 1))
  expect_named(an$y, c("y1", "y2"))
  expect_equal(colnames(an$sensitivity), c("a", "b", "c"))
})

testthat::test_that("the sensitivity is taken for the numeric elements of parms, or those asked for", {
  mixed <- list(a = 2.0, b = 1.5, c = 0.7, label = "demog", extra = c(1, 2))
  eq <- ode_steady_state(demog_rhs, c(0, 0), mixed)
  expect_equal(colnames(eq$sensitivity), c("a", "b", "c"))
  expect_equal(unname(eq$sensitivity), demog_sensitivity(demog_pars), tolerance = 1e-5)
  two <- ode_steady_state(demog_rhs, c(0, 0), mixed, sensitivity = c("c", "a"))
  expect_equal(colnames(two$sensitivity), c("c", "a"))
  expect_equal(unname(two$sensitivity), demog_sensitivity(demog_pars)[, c(3, 1)], tolerance = 1e-5)
  byi <- ode_steady_state(demog_rhs, c(0, 0), demog_pars, sensitivity = 2)
  expect_equal(colnames(byi$sensitivity), "b")
  none <- ode_steady_state(demog_rhs, c(0, 0), demog_pars, sensitivity = FALSE)
  expect_null(none$sensitivity)
  expect_equal(ode_counts(none)[["n_rhs"]] + 6, ode_counts(eq)[["n_rhs"]])
  positional <- function(t, y, p) demog_rhs(t, y, c(a = p[[1]], b = p[[2]], c = p[[3]]))
  unnamed <- ode_steady_state(positional, c(0, 0), unname(demog_pars))
  expect_equal(colnames(unnamed$sensitivity), c("p1", "p2", "p3"))
  expect_equal(unname(unnamed$sensitivity), demog_sensitivity(demog_pars), tolerance = 1e-5)
  # No parameters at all: nothing to be sensitive to.
  nop <- ode_steady_state(function(t, y, p) 1 - y, 0)
  expect_null(nop$sensitivity)
  expect_equal(nop$y, c(y1 = 1), tolerance = 1e-10)
  expect_error(ode_steady_state(demog_rhs, c(0, 0), mixed, sensitivity = "label"),
               "not among the numeric parameters")
})

testthat::test_that("Newton lands on the trivial repelling root; warmup gets past it", {
  bare <- ode_steady_state(logistic_rhs, c(1e-3, 0), logistic_pars)
  expect_true(bare$converged)
  expect_false(bare$warmed)
  expect_equal(unname(bare$y), c(0, 0), tolerance = 1e-10)
  expect_false(bare$stable)
  expect_equal(bare$spectral_abscissa, 0.6, tolerance = 1e-5)

  warm <- ode_steady_state(logistic_rhs, c(1e-3, 0), logistic_pars, warmup = 50, autonomous = TRUE)
  expect_true(warm$converged)
  expect_true(warm$warmed)
  expect_true(warm$stable)
  expect_equal(unname(warm$y), logistic_equilibrium(logistic_pars), tolerance = 1e-9)
  expect_equal(unname(warm$sensitivity), logistic_sensitivity(logistic_pars), tolerance = 1e-5)
  expect_equal(warm$time_dependence, 0)
  # The warm-up's evaluations are in the count.
  expect_gt(ode_counts(warm)[["n_rhs"]], ode_counts(bare)[["n_rhs"]] + 50)
  expect_output(print(warm), "after a warm-up")

  # Explicit landing times do the same.
  times <- ode_steady_state(logistic_rhs, c(1e-3, 0), logistic_pars, warmup = c(0, 25, 50))
  expect_equal(times$y, warm$y, tolerance = 1e-9)
  # A warm-up is not run when Newton already found an attractor.
  direct <- ode_steady_state(logistic_rhs, c(5, 5), logistic_pars, warmup = 50)
  expect_false(direct$warmed)
  expect_equal(direct$y, warm$y, tolerance = 1e-9)
})

testthat::test_that("a solve that does not converge warns and withholds what it cannot know", {
  expect_warning(eq <- ode_steady_state(demog_rhs, c(50, 50), demog_pars, max_iter = 1),
                 "did not converge")
  expect_false(eq$converged)
  expect_equal(eq$iterations, 1)
  expect_null(eq$jacobian)
  expect_null(eq$eigenvalues)
  expect_null(eq$stable)
  expect_null(eq$sensitivity)
  expect_true(is.na(eq$time_dependence))
  expect_output(print(eq), "NOT converged")
})

testthat::test_that("a singular Jacobian at the root is reported as a bifurcation point", {
  # f = k y^2 has a double root at 0, where df/dy = 0: the implicit function
  # theorem has nothing to say. jacfunc gives the exact zero.
  expect_error(ode_steady_state(function(t, y, p) p[["k"]] * y^2, 0, c(k = 1),
                                jacfunc = function(t, y, p) matrix(2 * p[["k"]] * y)),
               "singular")
})

testthat::test_that("explicit time dependence is reported, and skipped when declared autonomous", {
  forced <- function(t, y, p) 1 + t - y
  eq <- ode_steady_state(forced, 0, t0 = 2)
  expect_equal(eq$y, c(y1 = 3), tolerance = 1e-9)
  expect_equal(eq$time_dependence, 1, tolerance = 1e-5)
  expect_equal(ode_steady_state(forced, 0, t0 = 2, autonomous = TRUE)$time_dependence, 0)
})

testthat::test_that("bad arguments are refused before anything runs", {
  expect_error(ode_steady_state("f", 1), "func must be a function")
  expect_error(ode_steady_state(demog_rhs, numeric(0)), "y0 must be")
  expect_error(ode_steady_state(demog_rhs, c(0, 0), tol = -1), "tol must be")
  expect_error(ode_steady_state(demog_rhs, c(0, 0), warmup = -5), "warmup must be")
  expect_error(ode_steady_state(demog_rhs, c(0, 0), warmup = c(1, 2)), "start at t0")
  expect_error(ode_steady_state(demog_rhs, c(0, 0), control = 1, warmup = 1), "OdeControl")
  expect_error(ode_steady_state(function(t, y, p) stop("boom"), 1), "boom")
  expect_error(ode_steady_state(function(t, y, p) c(1, 2), 1), "returned 2 rates")
})

testthat::test_that("ode_steady_state agrees with rootSolve::steady on the same function", {
  skip_if_not_installed("rootSolve")
  ref <- rootSolve::steady(y = c(0, 0), func = demog_rhs, parms = demog_pars, method = "stode")
  eq <- ode_steady_state(demog_rhs, c(0, 0), demog_pars)
  expect_equal(unname(eq$y), unname(ref$y), tolerance = 1e-6)
})
