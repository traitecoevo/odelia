# $fit() on the leaf thermal example: a System read through time-varying drivers,
# so the reverse sweep has to stand each recorded row at its own time for the
# forcing to be right. The gradient is refereed against a central difference of
# the loss, which shares no code with the sweep.

leaf_fit_setup <- function() {
  ensure_leaf_thermal_interfaces(rebuild = FALSE)

  # Air temperature varying through the day, so the drivers matter at every step.
  time_driver <- seq(0, 24, by = 0.25)
  t_air <- 30 + 6 * sin(2 * pi * (time_driver - 15) / 24)
  drivers <- Drivers$new()
  drivers$set_variable("temperature", time_driver, t_air)

  pars_true <- list(k_H = 0.8, g_tr_max = 2.0, m_tr = 0.6, T_tr_mid = 28.0)
  sys_true <- LeafThermalSystem$new(pars_true, drivers)
  sys_true$set_initial_state(20.0, 0.0)
  ctrl <- OdeControl$new()
  runner <- LeafThermalSolver$new(sys_true$ptr, ctrl$ptr, drivers$ptr)
  # set_initial_state moves the reset point, not the current state, so reset
  # before the reference run: $fit() replays from the reset point.
  runner$reset()
  runner$advance_adaptive(seq(0, 12, by = 0.5))
  times <- runner$times()
  hist <- runner$history()

  pars_guess <- list(k_H = 0.5, g_tr_max = 1.0, m_tr = 0.5, T_tr_mid = 30.0)
  sys_fit <- LeafThermalSystem$new(pars_guess, drivers)
  sys_fit$set_initial_state(20.0, 0.0)
  fitter <- LeafThermalSolver$new(sys_fit$ptr, ctrl$ptr, drivers$ptr)
  fitter$set_target(times, matrix(hist$T_LC, ncol = 1), match(hist$time, times))

  list(fitter = fitter, truth = unlist(pars_true), guess = unlist(pars_guess))
}

central_difference <- function(f, x, rel = 1e-6) {
  vapply(seq_along(x), function(i) {
    h <- rel * max(1, abs(x[[i]]))
    up <- x; dn <- x
    up[[i]] <- up[[i]] + h
    dn[[i]] <- dn[[i]] - h
    (f(up) - f(dn)) / (2 * h)
  }, numeric(1))
}

testthat::test_that("leaf thermal fit is zero at the parameters that made the target", {
  s <- leaf_fit_setup()
  res <- s$fitter$fit(params = unname(s$truth))
  expect_lt(res$loss, 1e-20)
  expect_equal(res$gradient, rep(0, 4), tolerance = 1e-8)
})

testthat::test_that("leaf thermal parameter gradient matches a central difference", {
  s <- leaf_fit_setup()
  p <- unname(s$guess)
  res <- s$fitter$fit(params = p)
  expect_length(res$gradient, 4)
  fd <- central_difference(function(q) s$fitter$fit(params = q)$loss, p)
  expect_equal(res$gradient, fd, tolerance = 1e-6)
  expect_true(all(abs(res$gradient) > 0))
})

testthat::test_that("leaf thermal initial-state gradient matches a central difference", {
  s <- leaf_fit_setup()
  p <- unname(s$guess)
  res <- s$fitter$fit(ic = 25.0, params = p)
  expect_length(res$gradient, 5)
  fd <- central_difference(function(y) s$fitter$fit(ic = y, params = p)$loss, 25.0)
  expect_equal(res$gradient[5], fd, tolerance = 1e-6)
  # The parameter half does not depend on whether the initial state was asked for.
  expect_equal(res$gradient[1:4],
               s$fitter$fit(ic = 25.0, params = p)$gradient[1:4])
})

testthat::test_that("fit refuses a call it cannot answer", {
  s <- leaf_fit_setup()
  expect_error(s$fitter$fit(), "at least one of 'ic' or 'params'")
  expect_error(s$fitter$fit(params = c(1, 2)), "one entry per parameter")
})

testthat::test_that("a fit does not read the freed R system it was built from", {
  # leaf_fit_setup() drops its R-side systems on return, so the solver's own copy
  # must not point into them; if it does, a collection frees what it reads.
  s <- leaf_fit_setup()
  p <- unname(s$guess)
  before <- s$fitter$fit(ic = 25.0, params = p)
  invisible(gc(full = TRUE))
  invisible(gc(full = TRUE))
  expect_identical(s$fitter$fit(ic = 25.0, params = p), before)
})
