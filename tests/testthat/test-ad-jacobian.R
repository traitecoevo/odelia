# compute_jacobian differentiates a multi-output functional by delegating the
# record-once/row-sweep to xad::computeJacobian. Check the whole m x n Jacobian of the
# final state w.r.t. the parameters against central finite differences, and that the
# Solver-owned tape survives being reused across calls.

# Final state of a forward (double) run over `times`; the finite-difference
# reference for the differentiated multi-output functional.
final_state_vec <- function(pars, y0, times, ctrl) {
  sys <- LorenzSystem$new(pars[1], pars[2], pars[3])
  sys$set_initial_state(y0, t0 = times[1])
  solver <- Lorenz_Solver$new(sys$ptr, ctrl$ptr)
  solver$advance_fixed(times)
  solver$state()
}

test_that("compute_jacobian matches finite differences", {
  testthat::skip_if(is_pkgload_dll(), "Skipping AD workflow in pkgload load_all sessions due unstable native-pointer lifecycle.")

  pars <- c(sigma = 10.0, R = 28.0, b = 8.0 / 3.0)
  y0 <- c(1, 1, 1)
  times <- seq(0, 1, length.out = 11)
  ctrl <- odelia:::OdeControl$new()

  sys <- LorenzSystem$new(pars[1], pars[2], pars[3])
  sys$set_initial_state(y0, t0 = times[1])
  solver <- Lorenz_Solver$new(sys$ptr, ctrl$ptr)  # double handle; active built internally

  res <- Solver_jacobian_final_state(solver$ptr, times, pars)

  # Values reproduce the forward solve.
  expect_equal(res$values, final_state_vec(pars, y0, times, ctrl), tolerance = 1e-8)

  # Jacobian (m outputs x n params) matches central finite differences, one
  # column per parameter.
  eps <- 1e-6
  fd <- vapply(seq_along(pars), function(j) {
    hi <- lo <- pars
    hi[j] <- hi[j] + eps
    lo[j] <- lo[j] - eps
    (final_state_vec(hi, y0, times, ctrl) -
       final_state_vec(lo, y0, times, ctrl)) / (2 * eps)
  }, numeric(length(y0)))

  expect_equal(dim(res$jacobian), dim(fd))
  expect_equal(res$jacobian, unname(fd), tolerance = 1e-5)
})

test_that("repeated jacobian calls on one double solver reproduce", {
  testthat::skip_if(is_pkgload_dll(), "Skipping AD workflow in pkgload load_all sessions due unstable native-pointer lifecycle.")

  pars <- c(sigma = 10.0, R = 28.0, b = 8.0 / 3.0)
  y0 <- c(1, 1, 1)
  times <- seq(0, 1, length.out = 11)
  ctrl <- odelia:::OdeControl$new()

  sys <- LorenzSystem$new(pars[1], pars[2], pars[3])
  sys$set_initial_state(y0, t0 = times[1])
  solver <- Lorenz_Solver$new(sys$ptr, ctrl$ptr)  # double handle; active built internally

  first <- Solver_jacobian_final_state(solver$ptr, times, pars)
  second <- Solver_jacobian_final_state(solver$ptr, times, pars)

  # Each call builds its active replay and tape internally; a second call on the same
  # double solver must reproduce the first.
  expect_equal(second$values, first$values)
  expect_equal(second$jacobian, first$jacobian)
})
