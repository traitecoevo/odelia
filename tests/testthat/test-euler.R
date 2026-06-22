
# Fixed-step forward Euler is the alternative to the adaptive RKCK solver:
# a single derivative evaluation per step on a user-supplied grid (the way
# fixed-step DGVMs integrate). These tests exercise the Lorenz solver.

derivs_lorenz <- function(y, pars) {
  c(pars[[1]] * (y[[2]] - y[[1]]),
    pars[[2]] * y[[1]] - y[[2]] - y[[1]] * y[[3]],
    -pars[[3]] * y[[3]] + y[[1]] * y[[2]])
}

testthat::test_that("advance_euler matches a hand-rolled forward Euler", {
  pars <- c(sigma = 10.0, R = 28.0, b = 8.0 / 3.0)
  y0 <- c(1, 1, 1)
  times <- seq(0, 2, by = 0.001)

  lz <- LorenzSystem$new(pars[[1]], pars[[2]], pars[[3]])
  lz$set_state(y0, 0)
  ctrl <- OdeControl$new()
  runner <- Lorenz_Solver$new(lz$ptr, ctrl$ptr)
  runner$advance_euler(times)

  ## Reference: explicit forward Euler in R over the same grid.
  y <- y0
  for (i in 2:length(times)) {
    h <- times[[i]] - times[[i - 1]]
    y <- y + h * derivs_lorenz(y, pars)
  }

  testthat::expect_equal(runner$state(), y, tolerance = 1e-12)

  out <- runner$history()
  testthat::expect_equal(nrow(out), length(times))
  testthat::expect_equal(out$time, times)
})

testthat::test_that("advance_euler does plain Euler, not the RKCK step", {
  pars <- c(10.0, 28.0, 8.0 / 3.0)
  y0 <- c(1, 1, 1)
  times <- seq(0, 1, by = 0.05) # coarse: methods must visibly diverge

  ctrl <- OdeControl$new()

  lz1 <- LorenzSystem$new(pars[[1]], pars[[2]], pars[[3]]); lz1$set_state(y0, 0)
  euler <- Lorenz_Solver$new(lz1$ptr, ctrl$ptr)
  euler$advance_euler(times)

  lz2 <- LorenzSystem$new(pars[[1]], pars[[2]], pars[[3]]); lz2$set_state(y0, 0)
  adaptive <- Lorenz_Solver$new(lz2$ptr, ctrl$ptr)
  adaptive$advance_adaptive(times)

  ## One explicit Euler update by hand.
  y <- y0
  for (i in 2:length(times)) {
    h <- times[[i]] - times[[i - 1]]
    y <- y + h * derivs_lorenz(y, pars)
  }

  testthat::expect_equal(euler$state(), y, tolerance = 1e-12)
  testthat::expect_false(isTRUE(all.equal(euler$state(), adaptive$state())))
})

testthat::test_that("forward Euler converges to the adaptive solution as the step shrinks", {
  pars <- c(10.0, 28.0, 8.0 / 3.0)
  y0 <- c(1, 1, 1)
  t_end <- 1.0

  ## Tight adaptive reference.
  ctrl_ref <- OdeControl$new()
  ctrl_ref$set_tol_rel(1e-10)
  ctrl_ref$set_tol_abs(1e-10)
  lz_ref <- LorenzSystem$new(pars[[1]], pars[[2]], pars[[3]]); lz_ref$set_state(y0, 0)
  ref_runner <- Lorenz_Solver$new(lz_ref$ptr, ctrl_ref$ptr)
  ref_runner$advance_adaptive(c(0, t_end))
  ref <- ref_runner$state()

  ctrl <- OdeControl$new()
  err <- vapply(c(0.01, 0.005, 0.0025), function(dt) {
    lz <- LorenzSystem$new(pars[[1]], pars[[2]], pars[[3]]); lz$set_state(y0, 0)
    runner <- Lorenz_Solver$new(lz$ptr, ctrl$ptr)
    runner$advance_euler(seq(0, t_end, by = dt))
    max(abs(runner$state() - ref))
  }, numeric(1))

  ## Monotone first-order convergence.
  testthat::expect_true(all(diff(err) < 0))
})

testthat::test_that("advance_euler validates its time grid", {
  lz <- LorenzSystem$new(10, 28, 8 / 3); lz$set_state(c(1, 1, 1), 0)
  runner <- Lorenz_Solver$new(lz$ptr, OdeControl$new()$ptr)
  testthat::expect_error(runner$advance_euler(numeric(0)),
                         "must be vector of at least length 1")
  testthat::expect_error(runner$advance_euler(c(1.0, 2.0)),
                         "First element in 'times' must be same as current time")
})
