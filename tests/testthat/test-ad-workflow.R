# Using a stateful function to memoise the AD avoids re-running the solver to
# give optim the gradient
memoise_fit <- function(runner) {
  cache_p <- NULL
  cache_res <- NULL

  fn <- function(p) {
    cache_res <<- runner$fit(params = p)
    cache_p <<- p
    cache_res$loss
  }

  gr <- function(p) {
    if (!identical(p, cache_p)) fn(p) # Defensive: populate cache if stale
    cache_res$gradient
  }

  list(obj = fn, grad = gr)
}

test_that("AD workflow optimizes Lorenz parameters", {
  testthat::skip_if(is_pkgload_dll(), "Skipping AD workflow in pkgload load_all sessions due unstable native-pointer lifecycle.")

  true_pars <- c(sigma = 10.0, R = 28.0, b = 8.0 / 3.0)
  times <- seq(0, 1.5, length.out = 100)

  # Generate reference trajectory
  lz <- LorenzSystem$new(sigma = true_pars[1], R = true_pars[2], b = true_pars[3])
  lz$set_initial_state(c(1, 1, 1), t0 = 0)
  ctrl <- odelia:::OdeControl$new()

  runner <- Lorenz_Solver$new(lz$ptr, ctrl$ptr)
  runner$advance_fixed(times)
  hist <- runner$history()

  # Create AD solver
  target_vals <- as.matrix(hist[, c("x", "y", "z")])
  ad_runner <- Lorenz_Solver$new(lz$ptr, ctrl$ptr, active = TRUE)
  ad_runner$set_target(times, target_vals, seq_along(times))

  # Set parameters to something different
  initial_guess <- c(sigma = 12.0, R = 30.0, b = 3.0)
  lz$set_params(initial_guess)

  expect_false(isTRUE(all.equal(lz$pars(), true_pars)))

  # Optimize to recover the true parameters
  helper <- memoise_fit(ad_runner)
  res <- optim(
    par = initial_guess,
    fn = helper$obj,
    gr = helper$grad,
    method = "L-BFGS-B",
    control = list(maxit = 100)
  )

  expect_true(is.numeric(res$value) && is.finite(res$value))

  expect_equal(res$par[1], true_pars[1])
  expect_equal(res$par[2], true_pars[2])
  expect_equal(res$par[3], true_pars[3])

  # Verify optimized parameters reproduce the target trajectory
  lz$set_params(res$par)
  runner$reset()
  runner$advance_adaptive(times)
  optimized_hist <- runner$history()

  expect_equal(nrow(optimized_hist), nrow(hist))
  expect_equal(optimized_hist$x, hist$x, tolerance = 1e-4)
  expect_equal(optimized_hist$y, hist$y, tolerance = 1e-4)
  expect_equal(optimized_hist$z, hist$z, tolerance = 1e-4)
})
