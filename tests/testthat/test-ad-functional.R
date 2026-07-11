# Any caller-supplied functional can be differentiated, not just a fixed loss. Here a
# second functional -- the summed final state after a fixed-grid advance, with no
# observations -- is checked against central finite differences of the forward solve.
# With the least-squares path in test-ad-workflow.R this pins the seam as
# functional-agnostic.

# Summed final state of a forward (double) run over `times`; the
# finite-difference reference for the differentiated functional.
final_state_sum <- function(pars, y0, times, ctrl) {
  sys <- LorenzSystem$new(pars[1], pars[2], pars[3])
  sys$set_initial_state(y0, t0 = times[1])
  solver <- Lorenz_Solver$new(sys$ptr, ctrl$ptr)
  solver$advance_fixed(times)
  sum(solver$state())
}

test_that("compute_gradient differentiates a custom functional", {
  testthat::skip_if(is_pkgload_dll(), "Skipping AD workflow in pkgload load_all sessions due unstable native-pointer lifecycle.")

  pars <- c(sigma = 10.0, R = 28.0, b = 8.0 / 3.0)
  y0 <- c(1, 1, 1)
  times <- seq(0, 1, length.out = 11)
  ctrl <- odelia:::OdeControl$new()

  sys <- LorenzSystem$new(pars[1], pars[2], pars[3])
  sys$set_initial_state(y0, t0 = times[1])
  solver <- Lorenz_Solver$new(sys$ptr, ctrl$ptr)  # double handle; active built internally

  res <- Solver_gradient_final_state(solver$ptr, times, pars)

  # Value reproduces the forward solve.
  expect_equal(res$value, final_state_sum(pars, y0, times, ctrl), tolerance = 1e-8)

  # Gradient matches central finite differences (AD is exact; FD is the check).
  eps <- 1e-6
  fd <- vapply(seq_along(pars), function(i) {
    hi <- lo <- pars
    hi[i] <- hi[i] + eps
    lo[i] <- lo[i] - eps
    (final_state_sum(hi, y0, times, ctrl) -
       final_state_sum(lo, y0, times, ctrl)) / (2 * eps)
  }, numeric(1))

  expect_equal(res$gradient, unname(fd), tolerance = 1e-5)
})
