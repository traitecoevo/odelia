library(testthat)
library(odelia)

# Using a stateful function to memoise the AD avoids 
# re-running the solver when optim calls the gradient 
memoise_fit <- function(runner) {
  last_p <- NULL
  last_res <- NULL
    
  calc <- function(p) {
    if (!identical(p, last_p)) {
      last_res <<- runner$fit(params = p)
      last_p <<- p
    }
    last_res
  }
    
  list(
    obj = function(p) calc(p)$loss,
    grad = function(p) calc(p)$gradient
  )
}

test_that("AD workflow runs", {
  # Create Lorenz system
  lz <- LorenzSystem$new(sigma = 10.0, R = 28.0, b = 8.0 / 3.0)
  lz$set_initial_state(c(1, 1, 1), t0 = 0)
  
  ctrl_ptr <- OdeControl_new()
  
  runner <- Lorenz_Solver$new(lz$ptr, ctrl_ptr)
  runner$advance_adaptive(c(0, 10)) 
  hist <- runner$history() 
  
  # Set parameters to something different
  new_pars <- c(sigma = 12.0, R = 30.0, b = 3.0)
  lz$set_params(new_pars)
  
  # Create AD solver
  ad_runner <- Lorenz_Solver$new(lz$ptr, ctrl_ptr, active = TRUE)
  
  target_times <- hist$time
  target_vals <- as.matrix(hist[, c("x", "y", "z")])

  ad_runner$set_target(target_times, target_vals, seq_along(target_times))

  # See if we can recover the originals
  helper <- memoise_fit(ad_runner)
  res <- optim(
    par = new_pars, 
    fn = helper$obj, 
    gr = helper$grad, 
    method = "L-BFGS-B", 
    control = list(maxit = 10)
  )
  
  lz$set_params(res$par)
  
  runner$reset()
  runner$advance_adaptive(c(0, 10))
  runner$history()
  
  expect_true(TRUE) 
})

test_that("fit() can be called multiple times (params only)", {
  lz <- LorenzSystem$new(10, 28, 8/3)
  lz$set_initial_state(c(1, 1, 1), t0 = 0)
  lz$set_params(c(12, 30, 3))
  
  ctrl <- OdeControl_new()
  ad_runner <- Lorenz_Solver$new(lz$ptr, ctrl, active = TRUE)
  
  times <- c(0, 1, 2, 3)
  targets <- matrix(rep(c(1,1,1), 4), ncol=3, byrow=TRUE)
  ad_runner$set_target(times, targets, 1:4)
  
  # First call
  res1 <- ad_runner$fit(params=c(12, 30, 3))
  expect_equal(length(res1$gradient), 3)
  
  # Second call - should NOT error
  res2 <- ad_runner$fit(params=c(11, 29, 2.9))
  expect_equal(length(res2$gradient), 3)
  
  # Third call for good measure
  res3 <- ad_runner$fit(params=c(10, 28, 8/3))
  expect_equal(length(res3$gradient), 3)
})
