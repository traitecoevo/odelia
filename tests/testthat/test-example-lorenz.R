
testthat::test_that("lorenz system runs and produces expected results", {
  testthat::skip_on_cran()

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
testthat::test_that("lorenz AD IC gradients are NON-ZERO", {
  testthat::skip_on_cran()

  lz_true <- LorenzSystem$new(10, 28, 8/3)
  lz_true$set_initial_state(c(1, 1, 1), 0)

  ctrl <- OdeControl$new()
  runner_true <- Lorenz_Solver$new(lz_true$ptr, ctrl$ptr, active = FALSE)
  runner_true$advance_adaptive(seq(0, 5, by=0.25))
  hist_true <- runner_true$history()

  target_times <- hist_true$time
  target_vals <- as.matrix(hist_true[, c("x", "y", "z")])

  lz_fit <- LorenzSystem$new(10, 28, 8/3)
  lz_fit$set_initial_state(c(1, 1, 1), 0)

  ad_runner <- Lorenz_Solver$new(lz_fit$ptr, ctrl$ptr, active = TRUE)
  ad_runner$set_target(target_times, target_vals, c(1L, 2L, 3L))

  wrong_ic <- c(2, 2, 2)
  result <- ad_runner$fit(ic = wrong_ic, params = NULL)

  expect_true(is.finite(result$loss))
  expect_true(all(is.finite(result$gradient)))
  expect_equal(length(result$gradient), 3)
  
  expect_true(any(result$gradient != 0),
    info = "IC gradients are all zero - AD tape not capturing IC dependencies!")
})
testthat::test_that("lorenz AD parameter gradients are NON-ZERO", {
  testthat::skip_on_cran()

  true_pars <- c(sigma = 10.0, R = 28.0, b = 8.0/3.0)
  
  lz_true <- LorenzSystem$new(true_pars[1], true_pars[2], true_pars[3])
  lz_true$set_initial_state(c(1, 1, 1), 0)

  ctrl <- OdeControl$new()
  runner_true <- Lorenz_Solver$new(lz_true$ptr, ctrl$ptr, active = FALSE)
  runner_true$advance_adaptive(seq(0, 10, by=0.5))
  hist_true <- runner_true$history()

  target_times <- hist_true$time
  target_vals <- as.matrix(hist_true[, c("x", "y", "z")])

  wrong_pars <- c(sigma = 5.0, R = 20.0, b = 2.0)
  lz_fit <- LorenzSystem$new(wrong_pars[1], wrong_pars[2], wrong_pars[3])
  lz_fit$set_initial_state(c(1, 1, 1), 0)
  lz_fit$set_params(wrong_pars)

  ad_runner <- Lorenz_Solver$new(lz_fit$ptr, ctrl$ptr, active = TRUE)
  ad_runner$set_target(target_times, target_vals, c(1L, 2L, 3L))

  result <- ad_runner$fit(ic = NULL, params = wrong_pars)

  expect_true(is.finite(result$loss))
  expect_true(all(is.finite(result$gradient)))
  expect_equal(length(result$gradient), 3)
  
  expect_true(any(result$gradient != 0),
    info = "Lorenz parameter gradients should be non-zero (this is the only working AD gradient!)")
})