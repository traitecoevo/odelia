library(odelia)

setup_leaf_thermal_interfaces <- function(rebuild = FALSE) {
  withr::local_dir(here::here("inst/examples/leaf_thermal"))

  pkg_include <- here::here("inst/include")
  withr::local_envvar(PKG_CPPFLAGS = paste0("-I", pkg_include))

  expect_true(isNamespaceLoaded("odelia"))

  odelia_so <- file.path(
    system.file("libs", .Platform$r_arch, package = "odelia"),
    paste0("odelia", .Platform$dynlib.ext)
  )
  expect_true(file.exists(odelia_so))
  expect_silent(suppressWarnings(try(dyn.load(odelia_so, local = FALSE, now = TRUE), silent = TRUE)))

  expect_no_error(
    Rcpp::sourceCpp("../../../src/ode_interface.cpp", rebuild = rebuild, verbose = FALSE)
  )
  expect_no_error(
    Rcpp::sourceCpp("src/leaf_thermal_interface.cpp", rebuild = rebuild, verbose = FALSE)
  )
  expect_silent(source("R/leaf_thermal_interface.R"))
}

testthat::test_that("leaf thermal AD setup compiles", {
  testthat::skip_on_cran()
  setup_leaf_thermal_interfaces(rebuild = FALSE)
})

testthat::test_that("leaf thermal set_initial_state and reset work", {
  testthat::skip_on_cran()
  setup_leaf_thermal_interfaces(rebuild = FALSE)

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
  setup_leaf_thermal_interfaces(rebuild = FALSE)

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
  setup_leaf_thermal_interfaces(rebuild = FALSE)

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
  setup_leaf_thermal_interfaces(rebuild = FALSE)

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
