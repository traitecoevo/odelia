#!/usr/bin/env Rscript
# Test script for leaf thermal XAD integration

library(testthat)
library(odelia)

cat("Testing Leaf Thermal AD Integration\n")

# Set working directory to leaf thermal example folder
setwd(here::here("inst/examples/leaf_thermal"))

# 1. Test compilation
cat("1. Testing compilation...\n")
pkg_include <- here::here("inst/include")
Sys.setenv(PKG_CPPFLAGS = paste0("-I", pkg_include))

test_that("Leaf thermal compiles", {
  expect_silent(
    Rcpp::sourceCpp("src/leaf_thermal_interface.cpp", rebuild = TRUE, verbose = TRUE)
  )
})
cat("   ✓ Compilation successful\n\n")

# 2. Source R interface
cat("2. Loading R interface...\n")
source("R/leaf_thermal_interface.R")
cat("   ✓ R interface loaded\n\n")

# 3. Test basic (non-AD) simulation
cat("3. Testing basic simulation (double mode)...\n")
test_that("Basic leaf thermal simulation works", {
  # Create drivers
  drv <- Drivers$new()
  drv$set_constant("temperature", 25.0)
  
  # Create system with default parameters
  pars <- LeafThermalSystemPars()
  sys <- LeafThermalSystem$new(pars, drv)
  
  # Check system was created
  expect_false(is.null(sys$ptr))
  
  # Set initial state
  sys$set_state(20.0, 0.0)
  expect_equal(sys$state(), 20.0)
  
  # Create control and solver
  ctrl <- OdeControl$new()
  runner <- LeafThermalSolver$new(sys$ptr, ctrl$ptr, drv$ptr)
  
  # Run simulation
  times <- seq(0, 10, by = 0.1)
  runner$advance_adaptive(times)
  hist <- runner$history()
  
  # Check output
  expect_equal(nrow(hist), length(times))
  expect_true(all(c("time", "T_LC", "T_air", "dT_LC", "S_tr") %in% names(hist)))
  expect_true(all(is.finite(hist$T_LC)))
})
cat("   ✓ Basic simulation works\n\n")

# 4. Test set_initial_state and reset
cat("4. Testing set_initial_state and reset...\n")
test_that("set_initial_state and reset work", {
  drv <- Drivers$new()
  drv$set_constant("temperature", 25.0)
  
  pars <- LeafThermalSystemPars()
  sys <- LeafThermalSystem$new(pars, drv)
  
  # Set initial state
  LeafThermalSystem_set_initial_state(sys$ptr, 22.0, 0.0)
  
  # Reset should go back to initial state
  LeafThermalSystem_reset(sys$ptr)
  
  # State should be at initial
  state <- LeafThermalSystem_state(sys$ptr)
  expect_equal(state, 22.0)
})
cat("   ✓ set_initial_state and reset work\n\n")

# 5. Test AD solver creation
cat("5. Testing AD solver creation...\n")
test_that("AD solver can be created", {
  drv <- Drivers$new()
  drv$set_constant("temperature", 25.0)
  
  pars <- LeafThermalSystemPars()
  sys <- LeafThermalSystem$new(pars, drv)
  
  LeafThermalSystem_set_initial_state(sys$ptr, 20.0, 0.0)
  
  ctrl <- OdeControl$new()
  
  # Create AD solver with active=TRUE
  solver_ad <- Solver_new(sys$ptr, ctrl$ptr, drv$ptr, active = TRUE)
  
  expect_false(is.null(solver_ad))
})
cat("   ✓ AD solver creation works\n\n")

# 6. Test set_target
cat("6. Testing set_target...\n")
test_that("set_target works for AD solver", {
  drv <- Drivers$new()
  drv$set_constant("temperature", 25.0)
  
  pars <- LeafThermalSystemPars()
  sys <- LeafThermalSystem$new(pars, drv)
  
  LeafThermalSystem_set_initial_state(sys$ptr, 20.0, 0.0)
  
  ctrl <- OdeControl$new()
  solver_ad <- Solver_new(sys$ptr, ctrl$ptr, drv$ptr, active = TRUE)
  
  # Create target data (simple constant temperature)
  times <- c(0.0, 1.0, 2.0)
  target <- matrix(c(20.0, 21.0, 22.0), nrow = 3, ncol = 1)
  obs_indices <- as.integer(1)
  
  # Set target
  expect_silent(
    Solver_set_target(solver_ad, times, target, obs_indices, active = TRUE)
  )
})
cat("   ✓ set_target works\n\n")

# 7. Test Solver_fit (gradient computation)
cat("7. Testing Solver_fit (gradient computation)...\n")
test_that("Solver_fit computes gradients", {
  drv <- Drivers$new()
  drv$set_constant("temperature", 25.0)
  
  pars <- list(k_H = 0.5, g_tr_max = 1.0, m_tr = 0.5, T_tr_mid = 30.0)
  sys <- LeafThermalSystem$new(pars, drv)
  
  LeafThermalSystem_set_initial_state(sys$ptr, 20.0, 0.0)
  
  ctrl <- OdeControl$new()
  solver_ad <- Solver_new(sys$ptr, ctrl$ptr, drv$ptr, active = TRUE)
  
  # Create target
  times <- c(0.0, 1.0, 2.0, 3.0)
  target <- matrix(c(20.0, 21.0, 22.0, 23.0), nrow = 4, ncol = 1)
  obs_indices <- as.integer(1)
  
  Solver_set_target(solver_ad, times, target, obs_indices, active = TRUE)
  
  # Compute gradient w.r.t. initial condition
  result_ic <- Solver_fit(solver_ad, ic = 20.0, params = NULL)
  
  expect_true("loss" %in% names(result_ic))
  expect_true("gradient" %in% names(result_ic))
  expect_true(is.numeric(result_ic$loss))
  expect_equal(length(result_ic$gradient), 1)  # 1 IC
  
  # Compute gradient w.r.t. parameters
  params_vec <- c(0.5, 1.0, 0.5, 30.0)  # k_H, g_tr_max, m_tr, T_tr_mid
  result_params <- Solver_fit(solver_ad, ic = NULL, params = params_vec)
  
  expect_true(is.numeric(result_params$loss))
  expect_equal(length(result_params$gradient), 4)  # 4 params
  
  # Check gradients are non-zero
  expect_true(any(abs(result_ic$gradient) > 1e-10))
  expect_true(any(abs(result_params$gradient) > 1e-10))
})
cat("   ✓ Solver_fit computes gradients\n\n")
cat("All tests passed! \n")