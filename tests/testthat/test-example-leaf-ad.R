testthat::test_that("leaf thermal AD setup compiles", {
  testthat::skip_on_cran()
  ensure_leaf_thermal_interfaces(rebuild = FALSE)
})

testthat::test_that("leaf thermal set_initial_state and reset work", {
  testthat::skip_on_cran()
  ensure_leaf_thermal_interfaces(rebuild = FALSE)

  drv <- Drivers$new()
  drv$set_constant("temperature", 25.0)

  pars <- LeafThermalSystemPars()
  sys <- LeafThermalSystem$new(pars, drv)

  LeafThermalSystem_set_initial_state(sys$ptr, 22.0, 0.0)
  LeafThermalSystem_reset(sys$ptr)

  state <- LeafThermalSystem_state(sys$ptr)
  expect_equal(state, 22.0)
})

testthat::test_that("leaf thermal AD solver can be created", {
  testthat::skip_on_cran()
  ensure_leaf_thermal_interfaces(rebuild = FALSE)

  drv <- Drivers$new()
  drv$set_constant("temperature", 25.0)

  pars <- LeafThermalSystemPars()
  sys <- LeafThermalSystem$new(pars, drv)
  LeafThermalSystem_set_initial_state(sys$ptr, 20.0, 0.0)

  ctrl <- OdeControl$new()
  solver_ad <- LeafSolver_new(sys$ptr, ctrl$ptr, drv$ptr, active = TRUE)

  expect_false(is.null(solver_ad))
})

testthat::test_that("leaf thermal AD set_target works", {
  testthat::skip_on_cran()
  ensure_leaf_thermal_interfaces(rebuild = FALSE)

  drv <- Drivers$new()
  drv$set_constant("temperature", 25.0)

  pars <- LeafThermalSystemPars()
  sys <- LeafThermalSystem$new(pars, drv)
  LeafThermalSystem_set_initial_state(sys$ptr, 20.0, 0.0)

  ctrl <- OdeControl$new()
  solver_ad <- LeafSolver_new(sys$ptr, ctrl$ptr, drv$ptr, active = TRUE)

  times <- c(0.0, 1.0, 2.0)
  target <- matrix(c(20.0, 21.0, 22.0), nrow = 3, ncol = 1)
  obs_indices <- as.integer(1)

  expect_silent(
    LeafSolver_set_target(solver_ad, times, target, obs_indices, active = TRUE)
  )
})

testthat::test_that("leaf thermal AD fit computes finite outputs", {
  testthat::skip_on_cran()
  ensure_leaf_thermal_interfaces(rebuild = FALSE)

  drv <- Drivers$new()
  drv$set_constant("temperature", 25.0)

  pars <- list(k_H = 0.5, g_tr_max = 1.0, m_tr = 0.5, T_tr_mid = 30.0)
  sys <- LeafThermalSystem$new(pars, drv)
  LeafThermalSystem_set_initial_state(sys$ptr, 20.0, 0.0)

  ctrl <- OdeControl$new()
  solver_ad <- LeafSolver_new(sys$ptr, ctrl$ptr, drv$ptr, active = TRUE)

  times <- c(0.0, 1.0, 2.0, 3.0)
  target <- matrix(c(20.0, 21.0, 22.0, 23.0), nrow = 4, ncol = 1)
  obs_indices <- as.integer(1)

  LeafSolver_set_target(solver_ad, times, target, obs_indices, active = TRUE)

  result_ic <- LeafSolver_fit(solver_ad, ic = 20.0, params = NULL)
  expect_true("loss" %in% names(result_ic))
  expect_true("gradient" %in% names(result_ic))
  expect_true(is.numeric(result_ic$loss))
  expect_equal(length(result_ic$gradient), 1)

  params_vec <- c(0.5, 1.0, 0.5, 30.0)
  result_params <- LeafSolver_fit(solver_ad, ic = NULL, params = params_vec)
  expect_true(is.numeric(result_params$loss))
  expect_equal(length(result_params$gradient), 4)

  expect_true(is.finite(result_ic$loss))
  expect_true(all(is.finite(result_ic$gradient)))
  expect_true(is.finite(result_params$loss))
  expect_true(all(is.finite(result_params$gradient)))
})
testthat::test_that("leaf thermal AD parameter gradients are NON-ZERO", {
  testthat::skip_on_cran()
  ensure_leaf_thermal_interfaces(rebuild = FALSE)

  drv_true <- Drivers$new()
  drv_true$set_constant("temperature", 30.0)

  pars_true <- list(k_H = 0.8, g_tr_max = 2.0, m_tr = 0.6, T_tr_mid = 28.0)
  sys_true <- LeafThermalSystem$new(pars_true, drv_true)
  LeafThermalSystem_set_initial_state(sys_true$ptr, 20.0, 0.0)

  ctrl_true <- OdeControl$new()
  solver_true <- LeafSolver_new(sys_true$ptr, ctrl_true$ptr, drv_true$ptr, active = FALSE)
  LeafSolver_advance_adaptive(solver_true, seq(0, 20, by=1), active = FALSE)
  hist_true <- LeafSolver_get_history(solver_true, active = FALSE)

  target_times <- hist_true[[1]]
  target_vals <- matrix(hist_true[[2]], ncol=1)

  drv_fit <- Drivers$new()
  drv_fit$set_constant("temperature", 30.0)

  pars_wrong <- list(k_H = 0.5, g_tr_max = 1.0, m_tr = 0.5, T_tr_mid = 30.0)
  sys_fit <- LeafThermalSystem$new(pars_wrong, drv_fit)
  LeafThermalSystem_set_initial_state(sys_fit$ptr, 20.0, 0.0)

  ctrl_fit <- OdeControl$new()
  solver_fit_ad <- LeafSolver_new(sys_fit$ptr, ctrl_fit$ptr, drv_fit$ptr, active = TRUE)
  LeafSolver_set_target(solver_fit_ad, target_times, target_vals, as.integer(1), active = TRUE)

  params_vec <- c(0.5, 1.0, 0.5, 30.0)
  result <- LeafSolver_fit(solver_fit_ad, ic = NULL, params = params_vec)

  expect_true(is.finite(result$loss))
  expect_true(all(is.finite(result$gradient)))
  expect_equal(length(result$gradient), 4)
  
  expect_true(any(result$gradient != 0), 
    info = "Parameter gradients are all zero - AD tape not capturing parameter dependencies!")
})

testthat::test_that("leaf thermal AD IC gradients are NON-ZERO", {
  testthat::skip_on_cran()
  ensure_leaf_thermal_interfaces(rebuild = FALSE)

  drv <- Drivers$new()
  drv$set_constant("temperature", 30.0)

  pars <- list(k_H = 0.8, g_tr_max = 2.0, m_tr = 0.6, T_tr_mid = 28.0)
  sys_true <- LeafThermalSystem$new(pars, drv)
  LeafThermalSystem_set_initial_state(sys_true$ptr, 20.0, 0.0)

  ctrl <- OdeControl$new()
  solver_true <- LeafSolver_new(sys_true$ptr, ctrl$ptr, drv$ptr, active = FALSE)
  LeafSolver_advance_adaptive(solver_true, seq(0, 10, by=0.5), active = FALSE)
  hist_true <- LeafSolver_get_history(solver_true, active = FALSE)

  target_times <- hist_true[[1]]
  target_vals <- matrix(hist_true[[2]], ncol=1)

  sys_fit <- LeafThermalSystem$new(pars, drv)
  LeafThermalSystem_set_initial_state(sys_fit$ptr, 20.0, 0.0)

  solver_fit_ad <- LeafSolver_new(sys_fit$ptr, ctrl$ptr, drv$ptr, active = TRUE)
  LeafSolver_set_target(solver_fit_ad, target_times, target_vals, as.integer(1), active = TRUE)

  wrong_ic <- 25.0
  result <- LeafSolver_fit(solver_fit_ad, ic = wrong_ic, params = NULL)

  expect_true(is.finite(result$loss))
  expect_true(is.finite(result$gradient))
  expect_equal(length(result$gradient), 1)
  
  expect_true(result$gradient != 0,
    info = "IC gradient is zero - AD tape not capturing IC dependencies!")
})