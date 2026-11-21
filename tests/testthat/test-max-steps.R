# OdeControl's max_steps: the most accepted steps one interval of an adaptive
# advance may take. Unlimited unless set, counted per interval between two
# requested times and per event search, and refused with an error that names
# the count and where the solve had reached.

vdp <- function(t, y, eps) {
  list(c(y[2], ((1 - y[1]^2) * y[2] - y[1]) / eps))
}

testthat::test_that("the default is no limit, and get_max_steps says so", {
  ctrl <- OdeControl$new()
  expect_identical(ctrl$get_max_steps(), Inf)
  expect_identical(ctrl$get_controls()$max_steps, Inf)
  ctrl$set_max_steps(50)
  expect_identical(ctrl$get_max_steps(), 50)
  expect_identical(ctrl$get_controls()$max_steps, 50)
  ctrl$set_max_steps(Inf)
  expect_identical(ctrl$get_max_steps(), Inf)
  expect_error(ctrl$set_max_steps(0), "at least 1")
  expect_error(ctrl$set_max_steps(2.5), "at least 1")
  expect_error(ctrl$set_max_steps(-3), "at least 1")
})

testthat::test_that("an interval that needs more steps than allowed fails, naming the count", {
  times <- seq(0, 2, by = 0.5)
  # Stiff Van der Pol: the explicit stepper takes thousands of steps an
  # interval, RODAS tens.
  free <- ode_solve(vdp, c(2, 0), times, parms = 1e-4, method = "dopri", autonomous = TRUE)
  expect_gt(ode_counts(free)[["n_steps"]], 4 * 200)
  expect_error(
    ode_solve(vdp, c(2, 0), times, parms = 1e-4, method = "dopri", autonomous = TRUE,
              max_steps = 200),
    "Took 200 steps, the most allowed \\(max_steps\\), from t = 0 without reaching t = 0.5")
  capped <- ode_solve(vdp, c(2, 0), times, parms = 1e-4, method = "rodas", autonomous = TRUE,
                      max_steps = 1000)
  expect_equal(capped, ode_solve(vdp, c(2, 0), times, parms = 1e-4, method = "rodas",
                                 autonomous = TRUE))
  # The cap is per interval: a limit the first interval fits under but a
  # later one does not names the later one. Under dense output the later
  # interval starts where the step that passed the first time ended.
  expect_error(
    ode_solve(vdp, c(2, 0), c(0, 0.01, 2), parms = 1e-4, method = "dopri",
              autonomous = TRUE, max_steps = 400),
    "Took 400 steps, the most allowed \\(max_steps\\), from t = 0.01[0-9]* without reaching t = 2")
})

testthat::test_that("the cap holds under dense output and a control given directly", {
  ctrl <- OdeControl$new()
  ctrl$set_tol_rel(1e-6)
  ctrl$set_tol_abs(1e-6)
  ctrl$set_max_steps(200)
  times <- seq(0, 2, by = 0.5)
  expect_error(
    ode_solve(vdp, c(2, 0), times, parms = 1e-4, method = "dopri", autonomous = TRUE,
              control = ctrl, dense = TRUE),
    "the most allowed \\(max_steps\\)")
  # ode_solve's own max_steps is ignored with a control given.
  expect_error(
    ode_solve(vdp, c(2, 0), times, parms = 1e-4, method = "dopri", autonomous = TRUE,
              control = ctrl, max_steps = Inf),
    "the most allowed \\(max_steps\\)")
})

testthat::test_that("an event search counts against the cap too", {
  s <- OdeSolver$new(function(t, y) y, y0 = 1, autonomous = TRUE,
                     control = OdeControl$new()$set_max_steps(3))
  expect_error(s$advance_until(function(t, y) y[[1]] - 1e6, time_max = 20),
               "Took 3 steps, the most allowed \\(max_steps\\)")
  s <- OdeSolver$new(function(t, y) y, y0 = 1, autonomous = TRUE)
  expect_equal(s$advance_until(function(t, y) y[[1]] - 2, time_max = 20), log(2),
               tolerance = 1e-6)
})

testthat::test_that("a compiled solver honours the cap", {
  lz <- LorenzSystem$new(sigma = 10, R = 28, b = 8 / 3)
  lz$set_state(c(1, 1, 1))
  ctrl <- OdeControl$new()
  ctrl$set_max_steps(10)
  runner <- Lorenz_Solver$new(lz$ptr, ctrl$ptr)
  expect_error(runner$advance_adaptive(c(0, 10)),
               "Took 10 steps, the most allowed \\(max_steps\\), from t = 0 without reaching t = 10")
  ctrl$set_max_steps(Inf)
  runner <- Lorenz_Solver$new(lz$ptr, ctrl$ptr)
  expect_silent(runner$advance_adaptive(c(0, 10)))
})

testthat::test_that("a count a double cannot hold exactly is refused", {
  ctrl <- OdeControl$new()
  expect_error(ctrl$set_max_steps(1e300), "at most 2\\^53")
  expect_error(ctrl$set_max_steps(2^53 + 2), "at most 2\\^53")
  ctrl$set_max_steps(2^53)
  expect_identical(ctrl$get_max_steps(), 2^53)
})
