test_that("Extrinsic Drivers", {

  testthat::skip_on_cran()
  ensure_leaf_thermal_interfaces(rebuild = FALSE)
  
  # create Drivers object
  expect_silent(drv <- Drivers$new())

  # check methods
  expect_contains(names(drv), c("clone", "evaluate", "evaluate_range", "get_names", "initialize", "set_constant", "set_extrapolate", "set_variable"))

  # set constant and evaluate
  expect_silent(drv$set_constant("test", 1000))
  expect_equal(drv$evaluate("test", 1), 1000)
  expect_equal(drv$evaluate("test", 10), 1000)

  # set new value
  expect_silent(drv$set_constant("test", 1.5))
  expect_equal(drv$evaluate("test", 1), 1.5)
  expect_equal(drv$evaluate("test", 10), 1.5)

  # set a bunch of values
  for (i in 1:10) {
    expect_silent(drv$set_constant(letters[i], i))
  }
  for (i in 1:10) {
    expect_equal(drv$evaluate(letters[i], 1), i)
  }

  # set new values
  for (i in 1:10) {
    expect_silent(drv$set_constant(letters[i], i - 1))
  }
  for (i in 1:10) {
    expect_equal(drv$evaluate(letters[i], 1), i - 1)
  }

  # expected erors
  expect_error(drv$evaluate("wrong", 1))

  # Now set variable data and evaluate
  x <- 1:10
  y <- sin(x)
  expect_silent(drv$set_variable("sin", x, y))
  for (i in 1:10) {
    expect_equal(drv$evaluate("sin", i), y[i])
  }

  # extrapolation should error by default
  expect_silent(drv$evaluate("sin", 1))
  expect_error(drv$evaluate("sin", 0)) 
  expect_error(drv$evaluate("sin", 11))

  # enable extrapolation
  expect_silent(drv$set_extrapolate("sin", TRUE))
  expect_silent(drv$evaluate("sin", 0))
  expect_silent(drv$evaluate("sin", 11))

  # Splines require sensible data
  expect_silent(drv <- Drivers$new())
  # lengths of x and y must be equal
  expect_error(drv$set_variable("bad", c(1, 2, 3), c(1, 2))) 
  #insufficient number of points")
  expect_error(drv$set_variable("bad", numeric(0), numeric(0)))
  # ascending order
  expect_error(drv$set_variable("bad", c(3, 2, 4), c(1, 1, 1)))
  # unique points
  expect_error(drv$set_variable("bad", c(1, 1, 2), c(1, 1, 1)))
  
  # Spline contains correct data
  xx <- seq(0, 10, length.out = 100)
  yy <- sin(xx)
  expect_silent(drv$set_variable("sine", xx, yy))
  
  xx_cmp <- seq(0, 10, length.out = 42)
  # TODO:: enable access to data for testing
  # expect_identical(s$xy, cbind(xx, yy, deparse.level = 0))
  # expect_identical(c(s$min, s$max), range(xx))

  # Splines are accurate enough
  yy_C <- drv$evaluate_range("sine", xx_cmp)
  expect_equal(yy_C, sin(xx_cmp), tolerance = 1e-6)
})
