# Event location through OdeSolver$advance_until. That the state left at
# an event is the one a run asked to land there holds, to the bit, is proven in
# R-free C++ (tests/standalone/event_location.cpp); these tests prove the R
# adapter: the time returned, the state left, NA, and refusals.

oscillator <- function(t, y) c(y[[2]], -y[[1]])

testthat::test_that("advance_until stops where event changes sign, under every method", {
  for (method in c("dopri", "rkck", "rodas")) {
    s <- OdeSolver$new(oscillator, y0 = c(1, 0), method = method, autonomous = TRUE)
    t <- s$advance_until(function(t, y) y[[1]], time_max = 10)
    expect_equal(t, pi / 2, tolerance = 1e-6)
    expect_identical(s$time(), t)
    expect_equal(s$state(), c(0, -1), tolerance = 1e-6)
  }
})

testthat::test_that("a falling and a rising crossing are both found, one after the other", {
  s <- OdeSolver$new(oscillator, y0 = c(1, 0), autonomous = TRUE)
  t1 <- s$advance_until(function(t, y) y[[1]] + 0.5, time_max = 10)
  t2 <- s$advance_until(function(t, y) y[[1]] + 0.5 - 1e-3, time_max = 10)
  expect_equal(t1, acos(-0.5), tolerance = 1e-8)
  # From just below -0.5 + 1e-3, cos t falls to -1 and comes back up through
  # it: a rising crossing.
  expect_gt(t2, t1)
  expect_equal(t2, 2 * pi - acos(-0.5 + 1e-3), tolerance = 1e-8)
})

testthat::test_that("a series of levels is found from where each search stopped", {
  s <- OdeSolver$new(function(t, y) y, y0 = 1, autonomous = TRUE)
  levels <- c(2, 5, 10, 50)
  t <- vapply(levels, function(L) {
    s$advance_until(function(t, y) y - L, time_max = 100)
  }, numeric(1))
  expect_equal(t, log(levels), tolerance = 1e-8)
})

testthat::test_that("with no event before time_max it returns NA and stops there", {
  s <- OdeSolver$new(function(t, y) y, y0 = 1, autonomous = TRUE)
  r <- OdeSolver$new(function(t, y) y, y0 = 1, autonomous = TRUE)
  expect_identical(s$advance_until(function(t, y) y - 100, time_max = 2), NA_real_)
  r$advance_adaptive(c(0, 2))
  expect_identical(s$time(), 2)
  expect_identical(s$state(), r$state())
  expect_identical(s$times(), r$times())
})

testthat::test_that("an event of zero at the start is an event there", {
  s <- OdeSolver$new(function(t, y) y, y0 = 1, autonomous = TRUE)
  expect_identical(s$advance_until(function(t, y) y - 1, time_max = 2), 0)
  expect_identical(s$times(), 0)
})

testthat::test_that("bad events and arguments are refused", {
  s <- OdeSolver$new(function(t, y) y, y0 = 1, autonomous = TRUE)
  expect_error(s$advance_until(function(t, y) c(y, y), time_max = 2),
               "event returned 2 values, expected one")
  expect_error(s$advance_until(function(t, y) NaN, time_max = 2),
               "g is not finite")
  expect_error(s$advance_until(function(t, y) y - 2, time_max = 2, tol = 0),
               "tol must be a single positive number")
  expect_error(s$advance_until(function(t, y) y - 2, time_max = -1),
               "before the current time")
  expect_error(s$advance_until(function(t, y) y - 2, time_max = Inf),
               "must be finite")
  # An error raised by the event leaves the solver between steps, usable.
  expect_error(s$advance_until(function(t, y) stop("no event here"), time_max = 2),
               "no event here")
  expect_equal(s$advance_until(function(t, y) y - 2, time_max = 2), log(2),
               tolerance = 1e-8)
})

testthat::test_that("an error while the crossing is located takes the step back", {
  a <- OdeSolver$new(oscillator, y0 = c(1, 0), autonomous = TRUE)
  b <- OdeSolver$new(oscillator, y0 = c(1, 0), autonomous = TRUE)
  negatives <- 0
  # The first negative value is the end of the step that crossed; the next
  # call is inside it.
  fails_inside <- function(t, y) {
    if (y[[1]] < 0) {
      negatives <<- negatives + 1
      if (negatives > 1) stop("an event that fails")
    }
    y[[1]]
  }
  expect_error(a$advance_until(fails_inside, time_max = 10), "an event that fails")
  expect_lt(a$time(), pi / 2)
  expect_gt(a$state()[[1]], 0)
  x <- function(t, y) y[[1]]
  expect_identical(a$advance_until(x, time_max = 10),
                   b$advance_until(x, time_max = 10))
  expect_identical(a$state(), b$state())
  expect_identical(a$times(), b$times())
})
