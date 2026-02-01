library(testthat)
library(odelia)

# Setup: Link against odelia and compile Lorenz example
so_name <- paste0("odelia", .Platform$dynlib.ext)
lib_path <- system.file("libs", package = "odelia")
lib_file <- list.files(lib_path, pattern = so_name, recursive = TRUE, full.names = TRUE)[1]

# Use full path to .so
Sys.setenv(PKG_LIBS = paste0(lib_file, " -Wl,-rpath,", dirname(lib_file)))

# Compile and load
Rcpp::sourceCpp("../../inst/examples/lorenz/src/lorenz_interface.cpp", verbose = TRUE)
source("../../inst/examples/lorenz/R/lorenz_interface.R")

test_that("AD workflow optimizes Lorenz parameters", {
  true_pars <- c(sigma = 10.0, R = 28.0, b = 8.0 / 3.0)
  
  # Generate reference trajectory
  lz <- LorenzSystem$new(sigma = true_pars[1], R = true_pars[2], b = true_pars[3])
  lz$set_initial_state(c(1, 1, 1), t0 = 0)
  ctrl_ptr <- OdeControl_new()
  
  runner <- Lorenz_Solver$new(lz$ptr, ctrl_ptr)
  runner$advance_adaptive(seq(0, 10, by = 0.01)) 
  hist <- runner$history()
    
    
    # Set parameters to something different
  initial_guess <- c(sigma = 12.0, R = 30.0, b = 3.0)
  
  # Create a fresh solver for each step to avoid XAD tape reuse issues
  run_ad <- function(p) {
    lz$set_params(p)
    ad_runner <- Lorenz_Solver$new(lz$ptr, ctrl_ptr, active = TRUE)
    
    target_times <- hist$time
    target_vals <- as.matrix(hist[, c("x", "y", "z")])
    ad_runner$set_target(target_times, target_vals, seq_along(target_times))
    res <- ad_runner$fit(params = p)
    # Explicitly remove runner to trigger GC if needed, though scope exit handles it
    rm(ad_runner) 
    res
  }

  # Using a stateful function to memoise the AD avoids re-running the solver to
  # give optim the gradient 
  memoise_fit <- function() {
    cache_p <- NULL
    cache_res <- NULL
    
    calculate <- function(p) {
      if (!identical(p, cache_p)) {
        cache_res <<- run_ad(p)
        cache_p <<- p
      }
    }
    
    list(
      obj = function(p) { calculate(p); cache_res$loss },
      grad = function(p) { calculate(p); cache_res$gradient }
    )
  }
  
  # Optimize to recover the true parameters
  helper <- memoise_fit()
  
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
  expect_equal(res$par[3], true_pars[3], tolerance = 0.1)
  
  # Verify optimized parameters reproduce the target trajectory
  lz$set_params(res$par)
  runner$reset()
  runner$advance_adaptive(seq(0, 10, by = 0.01))
  optimized_hist <- runner$history()
  
  # Ensure distinct times
  
  expect_equal(nrow(optimized_hist), nrow(hist))
  expect_equal(optimized_hist$x, hist$x, tolerance = 0.05)
  expect_equal(optimized_hist$y, hist$y, tolerance = 0.05)  
  expect_equal(optimized_hist$z, hist$z, tolerance = 0.05)
})
