test_that("Extrinsic Drivers", {
  
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

test_that("Non-uniform spline grid evaluates correctly", {
  # The grids above (1:10, seq(...)) are equidistant and hit the O(1) uniform
  # fast path. This case uses quadratically-spaced knots (dense near 0, sparse
  # near 10; gaps vary ~70x) so the non-uniform index branch is exercised. The
  # proportional index guess assumes uniform spacing, so on this grid it lands
  # many segments away from the true one, stress-testing the nudge correction.
  drv <- Drivers$new()
  x <- 10 * (seq(0, 1, length.out = 40))^2
  y <- sin(x)
  expect_silent(drv$set_variable("nu", x, y))

  # A cubic spline interpolates its control points exactly, so evaluating at the
  # knots must return y. This fails if the segment index is wrong at any knot,
  # including the first/last (extrapolation-edge) knots.
  for (i in seq_along(x)) {
    expect_equal(drv$evaluate("nu", x[i]), y[i])
  }

  # Between knots it should still approximate the underlying smooth function.
  x_cmp <- seq(0, 10, length.out = 60)
  expect_equal(drv$evaluate_range("nu", x_cmp), sin(x_cmp), tolerance = 1e-2)

  # Out-of-range queries on a non-uniform grid: error by default, then allowed
  # once extrapolation is enabled (left tail below min, right tail above max).
  expect_error(drv$evaluate("nu", -0.001))
  expect_error(drv$evaluate("nu", 10.001))
  expect_silent(drv$set_extrapolate("nu", TRUE))
  expect_silent(drv$evaluate("nu", -0.001))
  expect_silent(drv$evaluate("nu", 10.001))
})
