library(odelia)

test_that("Drivers can be instantiated", {

  testthat::skip_on_cran()
  
  # set working directory to example folder
  withr::local_dir(here::here("inst/examples/leaf_thermal"))
  pkg_include <- here::here("inst/include")
  withr::local_envvar(PKG_CPPFLAGS = paste0("-I", pkg_include))

  # check compilation
  expect_silent(
    Rcpp::sourceCpp("src/leaf_thermal_interface.cpp", rebuild = FALSE, verbose = FALSE))

  expect_silent(
    source("R/leaf_thermal_interface.R")
  )
  
  expected <- list(
    tol_abs = 1e-8,
    tol_rel = 1e-8,
    a_y = 1.0,
    a_dydt = 0.0,
    step_size_min = 1e-8,
    step_size_max = 10.0,
    step_size_initial = 1e-6
  )
  keys <- sort(names(expected))

  expect_silent(
  ctrl <- OdeControl$new()
  )

#  expect_identical(sort(names(ctrl)), keys)
#  expect_identical(unclass(ctrl)[keys], expected[keys])
})

