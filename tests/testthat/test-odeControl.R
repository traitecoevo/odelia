testthat::test_that("OdeControl can be instantiated", {

  testthat::skip_on_cran()
  
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

})

