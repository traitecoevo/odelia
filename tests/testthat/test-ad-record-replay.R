# Record on a double Solver, replay on an active Solver. CanopySystem is the example
# that opts into the Replayable hooks (Lorenz and leaf_thermal don't), so these tests
# check the record -> replay path end to end: that a fixed-node replay reproduces the
# adaptive run and that its gradient matches finite differences, on a fresh recording
# built per call.

test_that("replay reproduces the adaptive value and matches finite differences", {
  testthat::skip_if(is_pkgload_dll(), "AD workflow needs the installed DLL, not load_all.")

  ctrl <- odelia:::OdeControl$new()
  gain <- 2.0; y0 <- 1.0; Tmax <- 2.0

  res <- Canopy_record_replay_gradient(
    Canopy_new(gain, y0), ctrl$ptr, Tmax, reuse_light = FALSE)

  # The run genuinely takes several adaptive steps -- the schedule that gets recorded.
  expect_gt(res$n_steps, 1)

  # Replaying advance_fixed on the recorded schedule reproduces the adaptive final
  # value, with the light recomputed with the active scalar on the recorded nodes.
  adaptive <- Canopy_adaptive_final(Canopy_new(gain, y0), ctrl$ptr, Tmax)
  expect_equal(res$value, adaptive, tolerance = 1e-4)

  # Its gradient matches central finite differences of the fully adaptive run.
  eps <- 1e-6
  hi <- Canopy_adaptive_final(Canopy_new(gain + eps, y0), ctrl$ptr, Tmax)
  lo <- Canopy_adaptive_final(Canopy_new(gain - eps, y0), ctrl$ptr, Tmax)
  fd <- (hi - lo) / (2 * eps)
  expect_equal(res$gradient[[1]], fd, tolerance = 1e-4)
})

test_that("reusing recorded light reproduces the value and zeros its derivative", {
  testthat::skip_if(is_pkgload_dll(), "AD workflow needs the installed DLL, not load_all.")

  ctrl <- odelia:::OdeControl$new()
  gain <- 2.0; y0 <- 1.0; Tmax <- 2.0

  res <- Canopy_record_replay_gradient(
    Canopy_new(gain, y0), ctrl$ptr, Tmax, reuse_light = TRUE)

  # Reusing the recorded light values reproduces the trajectory cheaply.
  adaptive <- Canopy_adaptive_final(Canopy_new(gain, y0), ctrl$ptr, Tmax)
  expect_equal(res$value, adaptive, tolerance = 1e-3)

  # gain enters only through the light; with the light reused as constants its
  # derivative is zero by construction.
  expect_equal(res$gradient[[1]], 0, tolerance = 1e-8)
})

test_that("repeated record -> replay calls on one system reproduce", {
  testthat::skip_if(is_pkgload_dll(), "AD workflow needs the installed DLL, not load_all.")

  ctrl <- odelia:::OdeControl$new()
  sys <- Canopy_new(2.0, 1.0)

  a <- Canopy_record_replay_gradient(sys, ctrl$ptr, 2.0, reuse_light = FALSE)
  b <- Canopy_record_replay_gradient(sys, ctrl$ptr, 2.0, reuse_light = FALSE)

  expect_equal(b$value, a$value)
  expect_equal(b$gradient, a$gradient)
})
