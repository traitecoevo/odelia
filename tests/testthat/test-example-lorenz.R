
testthat::test_that("lorenz system runs and produces expected results", {

  # function to compute Lorenz derivatives in R for comparison
  derivs_lorenz <- function(y, pars) {
    sigma <- pars[[1]]
    R <- pars[[2]]
    b <- pars[[3]]
    c(
      sigma * (y[[2]] - y[[1]]),
      R * y[[1]] - y[[2]] - y[[1]] * y[[3]],
      -b * y[[3]] + y[[1]] * y[[2]]
    )
  }

  derivs_lorenz_d <- function(t, y, pars) {
    list(derivs_lorenz(y, pars))
  }
  
  # check compilation
  # comparison data
  pars <- c(sigma=10.0,
            R=28.0,
            b=8.0 / 3.0)
  t0 <- 0.0
  tt <- seq(t0, t0+2, by=0.001)
  y <- c(21, 21, 21)

  expect_silent({
    lz <- LorenzSystem$new(sigma = pars[[1]], R = pars[[2]], b = pars[[3]])

    # Set initial state 
    lz$set_state(c(1, 1, 1), 0.0)

    # Check current state and rates
    lz$state()
    lz$rates()

    # Set up ODE solver
    expect_silent(ctrl <- odelia:::OdeControl$new())

    # Advance with adaptive time stepping to time 10000
    runner <- Lorenz_Solver$new(lz$ptr, ctrl$ptr)
    times <- seq(0, 100, by = 0.1)
    runner$advance_adaptive(times)
    out <- runner$history()
  })

  # test parameters
  expect_equal(names(pars), c("sigma", "R", "b"))
  expect_equal(lz$pars(), pars |> unlist() |> unname())
  # check methods
  expect_contains(names(lz), c("clone", "pars", "ptr", "rates", "set_state", "state"))

  # check output
  expect_true(all(c("time", "x", "y", "z", "dxdt", "dydt", "dzdt") %in% names(out)))
  expect_equal(nrow(out), length(times))
  expect_equal(out$time, times)
  expect_false(any(out[1, ] == out[nrow(out), ]))
  expect_true(all(is.finite(out$x)))
  expect_lt(max(out$x), 20)
  expect_gt(min(out$x), -20)

  ## Check state setting
  expect_silent(lz$set_state(y))
  expect_silent(state <- lz$state())
  expect_equal(length(state), 3)
  expect_equal(state, y)
  expect_identical(lz$rates(), derivs_lorenz(lz$state(), pars))
  expect_null(lz$ode_time)

  ## Change the state 
  y2 <- runif(3)
  lz$set_state(y2)
  expect_equal(lz$state(), y2)
  expect_equal(lz$rates(), derivs_lorenz(y2, pars))

  ## Run it again
  expect_silent({
    times <- out$time
    runner$set_state(c(1, 1, 1), 0.0)
    runner$advance_adaptive(times)
    out2 <- runner$history()
  })
  expect_equal(2*nrow(out), nrow(out2)) # history accumualtes unless reset
  expect_equal(out, out2[1:1001,])
  expect_equal(out, out2[1002:2002, ])
  
  # Todo: reset history and re-run to check identical output

  ## Run with different tolerance, should be different results
  ctrl2 <- odelia:::OdeControl$new()
  ctrl2$set_tol_rel(1e-12)
  expect_silent({
    times <- out$time
    runner <- Lorenz_Solver$new(lz$ptr, ctrl2$ptr)
    runner$set_state(c(1, 1, 1), 0.0)
    runner$advance_adaptive(times)
    out2 <- runner$history()
  })

  expect_false(all(out[,-1] == out2[,-1]))

})

# $fit() against a reference run observed at eleven output times: the reverse sweep
# is taken one interval between observations at a time, so several observations
# are what exercise it. Refereed against a central difference of the loss.
lorenz_fit_setup <- function() {
  lz <- LorenzSystem$new(10, 28, 8 / 3)
  lz$set_initial_state(c(1, 1, 1), 0)
  ctrl <- odelia:::OdeControl$new()
  runner <- Lorenz_Solver$new(lz$ptr, ctrl$ptr)
  runner$advance_adaptive(seq(0, 1, by = 0.1))
  times <- runner$times()
  hist <- runner$history()
  fitter <- Lorenz_Solver$new(lz$ptr, ctrl$ptr)
  fitter$set_target(times, as.matrix(hist[, c("x", "y", "z")]),
                    match(hist$time, times))
  fitter
}

lorenz_central_difference <- function(f, x) {
  vapply(seq_along(x), function(i) {
    h <- 1e-6 * max(1, abs(x[[i]]))
    up <- x; dn <- x
    up[[i]] <- up[[i]] + h
    dn[[i]] <- dn[[i]] - h
    (f(up) - f(dn)) / (2 * h)
  }, numeric(1))
}

testthat::test_that("lorenz fit's initial-state gradient matches a central difference", {
  fitter <- lorenz_fit_setup()
  ic <- c(1.1, 0.9, 1.05)
  res <- fitter$fit(ic = ic)
  expect_length(res$gradient, 3)
  fd <- lorenz_central_difference(function(y) fitter$fit(ic = y)$loss, ic)
  expect_equal(res$gradient, fd, tolerance = 1e-6)
  expect_equal(fitter$fit(ic = c(1, 1, 1))$loss, 0, tolerance = 1e-20)
})

testthat::test_that("lorenz fit's parameter gradient matches a central difference", {
  fitter <- lorenz_fit_setup()
  p <- c(11, 27, 2.5)
  res <- fitter$fit(params = p)
  expect_length(res$gradient, 3)
  fd <- lorenz_central_difference(function(q) fitter$fit(params = q)$loss, p)
  expect_equal(res$gradient, fd, tolerance = 1e-6)
  # Both at once: parameters first, then the initial state.
  both <- fitter$fit(ic = c(1.1, 0.9, 1.05), params = p)
  expect_length(both$gradient, 6)
})
