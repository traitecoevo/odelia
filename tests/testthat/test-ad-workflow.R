library(testthat)
library(odelia)

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
    if (!identical(p, cache_p)) fn(p)  # Defensive: populate cache if stale
    cache_res$gradient
  }
  
  list(obj = fn, grad = gr)
}

test_that("AD workflow optimizes Lorenz parameters", {
  true_pars <- c(sigma = 10.0, R = 28.0, b = 8.0 / 3.0)
  
  # Generate reference trajectory
  lz <- LorenzSystem$new(sigma = true_pars[1], R = true_pars[2], b = true_pars[3])
  lz$set_initial_state(c(1, 1, 1), t0 = 0)
  ctrl_ptr <- OdeControl_new()
  
  runner <- Lorenz_Solver$new(lz$ptr, ctrl_ptr)
  runner$advance_adaptive(c(0, 10)) 
  hist <- runner$history() 
  
  # Set parameters to something different
  initial_guess <- c(sigma = 12.0, R = 30.0, b = 3.0)
  lz$set_params(initial_guess)
  
  # Create AD solver
  ad_runner <- Lorenz_Solver$new(lz$ptr, ctrl_ptr, active = TRUE)
  
  target_times <- hist$time
  target_vals <- as.matrix(hist[, c("x", "y", "z")])
  ad_runner$set_target(target_times, target_vals, seq_along(target_times))

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
  
  expect_equal(res$par[1], true_pars[1], tolerance = 0.1)
  expect_equal(res$par[2], true_pars[2], tolerance = 0.1)
  expect_equal(res$par[3], true_pars[3], tolerance = 0.1)
  
  # Verify optimized parameters reproduce the target trajectory
  lz$set_params(res$par)
  runner$reset()
  runner$advance_adaptive(c(0, 10))
  optimized_hist <- runner$history()
  
  expect_equal(nrow(optimized_hist), nrow(hist))
  expect_equal(optimized_hist$x, hist$x, tolerance = 0.05)
  expect_equal(optimized_hist$y, hist$y, tolerance = 0.05)  
  expect_equal(optimized_hist$z, hist$z, tolerance = 0.05)
})
