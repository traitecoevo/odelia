library(odelia)

testthat::test_that("leaf thermal example runs", {
  
  
  testthat::skip_on_cran()
  
  # set working directory to example folder
  withr::local_dir(here::here("inst/examples/leaf_thermal"))

  # check compilation
  # Add package include path
  pkg_include <- here::here("inst/include")
  withr::local_envvar(PKG_CPPFLAGS = paste0("-I", pkg_include))

  expect_silent(
    Rcpp::sourceCpp("src/leaf_thermal_interface.cpp", rebuild = FALSE, verbose = FALSE)
  )
  expect_silent( source("R/leaf_thermal_interface.R") )

  # drivers
  p <- list(Tmean = 32, Tamp = 6, tpeak = 15)
  time_driver <- seq(0, 48, by = 0.25)
  t_air <- p$Tmean + p$Tamp * sin(2 * pi * (time_driver - p$tpeak) / 24)
  
  expect_silent(drivers <- Drivers$new())
  expect_silent(drivers$set_variable("temperature", time_driver, t_air))

  # LeafThermalSystem run
  expect_silent({
    pars <- LeafThermalSystemPars()
    lz <- LeafThermalSystem$new(pars, drivers)
    lz$set_state(c(25), 0)   
    ctrl <- OdeControl$new()
    runner <- LeafThermalSolver$new(lz$ptr, ctrl$ptr)
    times <- seq(0, 48, by = 0.5)
    runner$advance_adaptive(times)
    out <- runner$history() |>
      dplyr::mutate(time = times)
  })

  # test parameters
  expect_equal(names(pars), c("k_H", "g_tr_max", "m_tr", "T_tr_mid"))
  expect_equal(lz$pars(), pars |> unlist() |> unname())
  # check methods
  expect_contains(names(lz), c("clone", "get_current_drivers", "initialize", "initialize_drivers", "pars", "ptr", "rates", "set_state", "state"))

  # check output
  expect_true(all(c( "time", "T_LC", "T_air", "dT_LC", "S_tr" ) %in% names(out)))
  expect_equal(nrow(out), length(times))
  expect_true(all(is.finite(out$T_LC)))
  expect_lt(max(out$T_LC), 40)
  expect_gt(min(out$T_LC), 20)
})

