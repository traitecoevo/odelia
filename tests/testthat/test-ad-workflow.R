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

  # Generate reference trajectory
  true_pars <- c(sigma = 10.0, R = 28.0, b = 8.0 / 3.0)
  lz <- LorenzSystem$new(sigma = true_pars[1], R = true_pars[2], b = true_pars[3])
  lz$set_initial_state(c(1, 1, 1), t0 = 0)
  ctrl <- odelia:::OdeControl$new()

  runner <- Lorenz_Solver$new(lz$ptr, ctrl$ptr)
  runner$advance_adaptive(c(0, 1))

  times <- runner$times()
  hist <- runner$history()
  obs_index <- which(times %in% hist$time)

  # Set parameters to something different
  initial_guess <- c(sigma = 12.0, R = 30.0, b = 3.0)
  lz$set_params(initial_guess)

  # A solver to fit with
  target_vals <- as.matrix(hist[, c("x", "y", "z")])
  ad_runner <- Lorenz_Solver$new(lz$ptr, ctrl$ptr)
  ad_runner$set_target(times, target_vals, obs_index)

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

  expect_equal(res$convergence, 0)
  expect_lt(res$value, 1e-10)

  expect_equal(res$par[1], true_pars[1], tolerance = 1e-6)
  expect_equal(res$par[2], true_pars[2], tolerance = 1e-6)
  expect_equal(res$par[3], true_pars[3], tolerance = 1e-6)

  # Verify optimized parameters reproduce the target trajectory
  lz$set_params(res$par)
  runner <- Lorenz_Solver$new(lz$ptr, ctrl$ptr)
  runner$advance_adaptive(c(0, 1))
  optimized_hist <- runner$history()

  expect_equal(nrow(optimized_hist), nrow(hist))
  expect_equal(optimized_hist$x, hist$x, tolerance = 1e-6)
  expect_equal(optimized_hist$y, hist$y, tolerance = 1e-6)
  expect_equal(optimized_hist$z, hist$z, tolerance = 1e-6)
})

test_that("a Dormand-Prince solver fits under Dormand-Prince", {
  testthat::skip_if(is_pkgload_dll(), "Skipping AD workflow in pkgload load_all sessions due unstable native-pointer lifecycle.")

  true_pars <- c(10.0, 28.0, 8.0 / 3.0)
  lz <- LorenzSystem$new(sigma = true_pars[1], R = true_pars[2], b = true_pars[3])
  lz$set_initial_state(c(1, 1, 1), t0 = 0)
  ctrl <- odelia:::OdeControl$new()
  reference <- Lorenz_Solver$new(lz$ptr, ctrl$ptr, method = "dopri")
  reference$advance_adaptive(c(0, 1))
  times <- reference$times()
  hist <- reference$history()
  obs_index <- which(times %in% hist$time)
  target <- as.matrix(hist[, c("x", "y", "z")])

  fitter <- function(method) {
    s <- Lorenz_Solver$new(lz$ptr, ctrl$ptr, method = method)
    s$set_target(times, target, obs_index)
    s
  }
  # At the parameters that made the target, the replay is the reference's own
  # integrator on its own steps, so the loss is rounding (about 1e-27); a
  # Cash-Karp replay of the same steps is a different integrator (about 3e-14).
  expect_lt(fitter("dopri")$fit(params = true_pars)$loss, 1e-20)
  expect_gt(fitter("rkck")$fit(params = true_pars)$loss, 1e-18)

  # And its gradient is the derivative of its own loss.
  s <- fitter("dopri")
  at <- c(10.5, 27.5, 2.8)
  eps <- 1e-6
  fd <- sapply(seq_along(at), function(j) {
    up <- at; up[j] <- up[j] + eps
    down <- at; down[j] <- down[j] - eps
    (s$fit(params = up)$loss - s$fit(params = down)$loss) / (2 * eps)
  })
  expect_equal(s$fit(params = at)$gradient, fd, tolerance = 1e-6)
})
