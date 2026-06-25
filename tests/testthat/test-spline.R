# Tests for Spline::deriv / Interpolator::deriv (the analytic first derivative).
# The interpolator/spline headers are header-only, so we just need the include
# path (no linking against the odelia shared library).

compile_spline_deriv_interface <- function() {
  include_dir <- dirname(dirname(resolve_test_path(
    "include/odelia/ode_solver.hpp", "inst/include/odelia/ode_solver.hpp")))
  withr::local_envvar(PKG_CPPFLAGS = paste0("-I", shQuote(include_dir)))
  Rcpp::sourceCpp(code = '
    #include <Rcpp.h>
    #include <vector>
    #include <odelia/interpolator.hpp>

    // [[Rcpp::export]]
    Rcpp::List spline_value_and_deriv(std::vector<double> x,
                                      std::vector<double> y,
                                      std::vector<double> u) {
      odelia::interpolator::Interpolator s;
      s.init(x, y);
      std::vector<double> value(u.size()), deriv(u.size());
      for (size_t i = 0; i < u.size(); ++i) {
        value[i] = s.eval(u[i]);
        deriv[i] = s.deriv(u[i]);
      }
      return Rcpp::List::create(Rcpp::_["value"] = value,
                                Rcpp::_["deriv"] = deriv);
    }', verbose = FALSE)
}

testthat::test_that("Interpolator::deriv reproduces a known cubic's derivative", {
  compile_spline_deriv_interface()

  # A cubic spline through samples of a cubic should recover the cubic exactly,
  # so its derivative should equal the analytic derivative 3x^2 + ... .
  f      <- function(x) 2 * x^3 - x^2 + 0.5 * x + 3
  fprime <- function(x) 6 * x^2 - 2 * x + 0.5
  x <- seq(-3, 3, length.out = 80)
  r <- spline_value_and_deriv(x, f(x), u <- c(-2.4, -0.7, 0.3, 1.9))

  expect_equal(r$value, f(u), tolerance = 1e-8)
  expect_equal(r$deriv, fprime(u), tolerance = 1e-6)
})

testthat::test_that("Interpolator::deriv agrees with a finite difference of eval", {
  compile_spline_deriv_interface()

  # For a smooth non-polynomial the spline is an approximation, but deriv() must
  # still be the exact derivative of whatever eval() interpolates -- so it agrees
  # with a central finite difference of eval() to high precision everywhere,
  # including across knots.
  x <- seq(0, 2 * pi, length.out = 200)
  y <- sin(x)
  u <- seq(0.3, 6.0, by = 0.37)
  r  <- spline_value_and_deriv(x, y, u)

  h  <- 1e-5
  rp <- spline_value_and_deriv(x, y, u + h)
  rm <- spline_value_and_deriv(x, y, u - h)
  fd <- (rp$value - rm$value) / (2 * h)

  expect_equal(r$deriv, fd, tolerance = 1e-5)
  # and it should track the true derivative cos(x) to spline accuracy
  expect_equal(r$deriv, cos(u), tolerance = 1e-3)
})
