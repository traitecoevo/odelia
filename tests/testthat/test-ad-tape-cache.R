# The active twin (and its tape) are cached on a persistent double Solver and reused
# across gradient calls, while the recording is read fresh on every call. Where
# test-ad-record-replay checks a single record -> replay is correct, these tests check
# the *reuse*: that reusing the cached twin reproduces the first result, that a new
# recording is picked up rather than replayed against a stale schedule, and that the
# two replay modes share one twin. CanopySystem is used because it is the example that
# records a schedule to go stale.

test_that("reused twin reproduces, and a fresh recording is picked up per call", {
  testthat::skip_if(is_pkgload_dll(), "AD workflow needs the installed DLL, not load_all.")

  ctrl <- odelia:::OdeControl$new()
  d <- Canopy_Solver_new(Canopy_new(2.0, 1.0), ctrl$ptr)

  # Record to Tmax = 2, then two gradient calls: the second reuses the cached twin
  # and tape and must reproduce the first exactly.
  n1 <- Canopy_record(d, Tmax = 2.0)
  expect_gt(n1, 1)
  g1  <- Canopy_replay_gradient(d, reuse_light = FALSE)
  g1b <- Canopy_replay_gradient(d, reuse_light = FALSE)
  expect_equal(g1b$value, g1$value)
  expect_equal(g1b$gradient, g1$gradient)

  # A fresh double pass (a longer schedule) is picked up on the next call -- the
  # reused twin is handed the new recording, not replayed against the stale Tmax = 2
  # one. The value matches a fresh adaptive run to 3.0.
  n2 <- Canopy_record(d, Tmax = 3.0)
  expect_gt(n2, n1)
  g2 <- Canopy_replay_gradient(d, reuse_light = FALSE)
  expect_false(isTRUE(all.equal(g2$value, g1$value)))

  adaptive3 <- Canopy_adaptive_final(Canopy_new(2.0, 1.0), ctrl$ptr, 3.0)
  expect_equal(g2$value, adaptive3, tolerance = 1e-4)
})

test_that("recompute and reuse-light replays share one cached twin, chosen per call", {
  testthat::skip_if(is_pkgload_dll(), "AD workflow needs the installed DLL, not load_all.")

  ctrl <- odelia:::OdeControl$new()
  d <- Canopy_Solver_new(Canopy_new(2.0, 1.0), ctrl$ptr)
  Canopy_record(d, Tmax = 2.0)

  # Recompute then reuse-light on the SAME cached twin: the mode is a per-call choice.
  recompute <- Canopy_replay_gradient(d, reuse_light = FALSE)
  reuse     <- Canopy_replay_gradient(d, reuse_light = TRUE)

  # Both reproduce the trajectory; reusing the light zeros its derivative, while
  # recomputing carries the real self-shading gradient.
  expect_equal(reuse$value, recompute$value, tolerance = 1e-3)
  expect_equal(reuse$gradient[[1]], 0, tolerance = 1e-8)
  expect_gt(abs(recompute$gradient[[1]]), 1e-6)
})

test_that("replay before record is a clear error, not a stale read", {
  testthat::skip_if(is_pkgload_dll(), "AD workflow needs the installed DLL, not load_all.")

  ctrl <- odelia:::OdeControl$new()
  d <- Canopy_Solver_new(Canopy_new(2.0, 1.0), ctrl$ptr)
  expect_error(Canopy_replay_gradient(d), "record")
})
