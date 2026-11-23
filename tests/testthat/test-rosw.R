# The Rosenbrock-W stepper (method = "ros34pw2") through R. Its order, and how
# many Jacobians it forms on stiff problems against RODAS, are tested in
# tests/standalone/rosw.cpp.

robertson_rhs <- function(t, y, p) {
  list(c(-0.04 * y[1] + 1e4 * y[2] * y[3],
         0.04 * y[1] - 1e4 * y[2] * y[3] - 3e7 * y[2]^2,
         3e7 * y[2]^2))
}
robertson_jac <- function(t, y, p) {
  matrix(c(-0.04, 1e4 * y[3], 1e4 * y[2],
           0.04, -1e4 * y[3] - 6e7 * y[2], -1e4 * y[2],
           0, 6e7 * y[2], 0), 3, 3, byrow = TRUE)
}

# A kept Jacobian that gets a step rejected is replaced, so the saving shrinks
# as the tolerance tightens: on this run RODAS forms 73 Jacobians at 1e-4 and
# ROS34PW2 29, against 187 and 131 at 1e-6.
testthat::test_that("ROS34PW2 solves Robertson's kinetics through R with fewer Jacobians than RODAS", {
  skip_if_not_installed("deSolve")
  times <- c(0, 10^seq(-2, 4, by = 0.5))
  run <- function(method) {
    ode_solve(robertson_rhs, c(1, 0, 0), times, jacfunc = robertson_jac,
              method = method, rtol = 1e-4, atol = 1e-10, autonomous = TRUE,
              dense = FALSE)
  }
  w <- run("ros34pw2")
  r <- run("rodas")
  ref <- deSolve::ode(c(1, 0, 0), times, robertson_rhs, NULL, method = "radau",
                      rtol = 1e-10, atol = 1e-14)
  expect_equal(w[, 2], as.numeric(ref[, 2]), tolerance = 1e-3)
  expect_equal(w[, 4], as.numeric(ref[, 4]), tolerance = 1e-3)
  expect_equal(rowSums(w[, 2:4]), rep(1, nrow(w)), tolerance = 1e-9)
  expect_lt(2 * ode_counts(w)[["n_jac"]], ode_counts(r)[["n_jac"]])
})

testthat::test_that("jacobian_max_age is how many steps one Jacobian serves", {
  ctrl <- OdeControl$new()
  expect_equal(ctrl$get_jacobian_max_age(), 20)
  ctrl$set_jacobian_max_age(Inf)
  expect_equal(ctrl$get_jacobian_max_age(), Inf)
  expect_equal(ctrl$get_controls()$jacobian_max_age, Inf)
  expect_error(ctrl$set_jacobian_max_age(0), "whole number of at least 1")
  expect_error(ctrl$set_jacobian_max_age(2.5), "whole number of at least 1")

  # At 1 a Jacobian is formed at every accepted step and never for a retry,
  # which starts where the rejected attempt did.
  ctrl$set_jacobian_max_age(1)
  ctrl$set_tol_rel(1e-6)
  ctrl$set_tol_abs(1e-12)
  s <- OdeSolver$new(robertson_rhs, c(1, 0, 0), jac = robertson_jac,
                     control = ctrl, method = "ros34pw2", autonomous = TRUE,
                     parms = list())
  s$advance_adaptive(c(0, 40))
  every <- s$counts()
  expect_equal(every[["n_jac"]], every[["n_steps"]])

  ctrl$set_jacobian_max_age(20)
  s <- OdeSolver$new(robertson_rhs, c(1, 0, 0), jac = robertson_jac,
                     control = ctrl, method = "ros34pw2", autonomous = TRUE,
                     parms = list())
  s$advance_adaptive(c(0, 40))
  kept <- s$counts()
  expect_lt(kept[["n_jac"]], every[["n_jac"]] / 2)
  expect_lte(kept[["n_jac"]], ceiling(kept[["n_steps"]] / 20) + kept[["n_rejections"]])
})

testthat::test_that("a solve under ROS34PW2 is not differentiated", {
  testthat::skip_if(is_pkgload_dll(), "Skipping AD workflow in pkgload load_all sessions due unstable native-pointer lifecycle.")
  lz <- LorenzSystem$new(sigma = 10, R = 28, b = 8 / 3)
  lz$set_initial_state(c(1, 1, 1), t0 = 0)
  ctrl <- OdeControl$new()
  runner <- Lorenz_Solver$new(lz$ptr, ctrl$ptr, method = "ros34pw2")
  expect_error(runner$sensitivities(c(0, 0.5)), "ros34pw2")
  runner$set_target(c(0, 0.25, 0.5), matrix(1, 2, 3), c(2L, 3L))
  expect_error(runner$loss_and_gradient(params = c(10, 28, 8 / 3)),
               "cannot be differentiated")
})
