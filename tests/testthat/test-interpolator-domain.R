# Tests for the out-of-domain message from Interpolator::eval when extrapolation
# is disabled. Header-only like test-spline.R, so only the include path is needed.

compile_interpolator_domain_interface <- function() {
  include_dir <- dirname(dirname(resolve_test_path(
    "include/odelia/ode_solver.hpp", "inst/include/odelia/ode_solver.hpp")))
  withr::local_envvar(PKG_CPPFLAGS = paste0("-I", shQuote(include_dir)))
  Rcpp::sourceCpp(code = '
    #include <Rcpp.h>
    #include <vector>
    #include <odelia/interpolator.hpp>

    // [[Rcpp::export]]
    double interpolator_eval_no_extrapolate(std::vector<double> x,
                                            std::vector<double> y,
                                            double u) {
      odelia::interpolator::Interpolator s;
      s.init(x, y);
      s.set_extrapolate(false);
      return s.eval(u);
    }', verbose = FALSE)
}

# A domain deliberately not [0, 1]: an endpoint of 1 would be indistinguishable
# from a hard-coded default in the message, and 6.8918 is the order of the
# hydraulic domain the message was written for (plant#576).
dom_x <- seq(0, 6.8918, length.out = 40)

testthat::test_that("an out-of-domain eval reports the point, the miss and the domain", {
  compile_interpolator_domain_interface()

  msg <- tryCatch(
    interpolator_eval_no_extrapolate(dom_x, dom_x^2, 7.5),
    error = conditionMessage)

  # Naming all three is the whole point: without them the message is identical
  # whichever spline threw it, and a consumer holding several learns nothing.
  expect_match(msg, "u = 7.5", fixed = TRUE)
  expect_match(msg, "beyond the upper end", fixed = TRUE)
  expect_match(msg, "0.6082", fixed = TRUE)        # how far out it fell
  expect_match(msg, "[0, 6.8918]", fixed = TRUE)   # the domain itself

  # Which end matters as much as the fact of missing: plant#576 was assumed to be
  # a point past the far end of a hydraulic curve, and was the near end. Aiming a
  # fix at the wrong end is the failure this half of the message prevents.
  msg_lo <- tryCatch(
    interpolator_eval_no_extrapolate(dom_x, dom_x^2, -0.0023),
    error = conditionMessage)
  expect_match(msg_lo, "beyond the lower end", fixed = TRUE)
  expect_match(msg_lo, "u = -0.0023", fixed = TRUE)
})

testthat::test_that("a non-finite eval point still falls through to the spline", {
  compile_interpolator_domain_interface()

  # Load-bearing, not an oversight. The guard is `u < min() || u > max()`, and
  # every comparison against NaN is false, so a non-finite u reaches the spline
  # and comes back non-finite. Consumers depend on it: plant documents a
  # profit_psi_stem_TF(NA, .) -> NA contract built on exactly this. Rewriting the
  # guard as `!(u >= min() && u <= max())` reads as a tightening and would turn
  # this into a throw.
  expect_true(is.nan(interpolator_eval_no_extrapolate(dom_x, dom_x^2, NaN)))
  expect_false(is.finite(interpolator_eval_no_extrapolate(dom_x, dom_x^2, NA_real_)))
})

testthat::test_that("an in-domain eval is unaffected", {
  compile_interpolator_domain_interface()

  expect_equal(interpolator_eval_no_extrapolate(dom_x, dom_x^2, 3.1),
               3.1^2, tolerance = 1e-6)
  # The endpoints are in the domain, not outside it: the comparison is strict.
  expect_equal(interpolator_eval_no_extrapolate(dom_x, dom_x^2, 0), 0,
               tolerance = 1e-9)
  expect_equal(interpolator_eval_no_extrapolate(dom_x, dom_x^2, 6.8918),
               6.8918^2, tolerance = 1e-6)
})
