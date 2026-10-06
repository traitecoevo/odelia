
testthat::test_that("leaf thermal example runs", {
  
  ensure_leaf_thermal_interfaces(rebuild = FALSE)

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
    runner <- LeafThermalSolver$new(lz$ptr, ctrl$ptr, drivers$ptr)
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

testthat::test_that("the leaf thermal example, a compiled system with drivers, steps under Dormand-Prince", {
  ensure_leaf_thermal_interfaces(rebuild = FALSE)
  p <- list(Tmean = 32, Tamp = 6, tpeak = 15)
  time_driver <- seq(0, 48, by = 0.25)
  drivers <- Drivers$new()
  drivers$set_variable("temperature", time_driver, p$Tmean + p$Tamp * sin(2 * pi * (time_driver - p$tpeak) / 24))
  times <- seq(0, 48, by = 0.5)
  run <- function(method) {
    lz <- LeafThermalSystem$new(LeafThermalSystemPars(), drivers)
    lz$set_state(c(25), 0)
    ctrl <- OdeControl$new()
    ctrl$set_tol_rel(1e-8)
    ctrl$set_tol_abs(1e-8)
    runner <- LeafThermalSolver$new(lz$ptr, ctrl$ptr, drivers$ptr, method = method)
    runner$advance_adaptive(times)
    runner$history()$T_LC
  }
  expect_equal(run("dopri"), run("rkck"), tolerance = 1e-6)
  # No Jacobian hook and no rebind() on this system, so the implicit stepper
  # says so rather than failing somewhere inside.
  lz <- LeafThermalSystem$new(LeafThermalSystemPars(), drivers)
  lz$set_state(c(25), 0)
  runner <- LeafThermalSolver$new(lz$ptr, OdeControl$new()$ptr, drivers$ptr, method = "rodas")
  expect_error(runner$advance_adaptive(c(0, 1)), "ode_jacobian\\(\\) hook")
})
