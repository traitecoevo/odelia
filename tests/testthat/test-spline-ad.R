# Tests for the scalar-templated spline (odelia::spline::basic_spline<S>), the
# differentiable-spline primitive for reverse-mode AD (#472 scope B /
# traitecoevo/plant#537). With knot POSITIONS x frozen (double) and knot VALUES
# active, the cubic-coefficient solve is a constant double band matrix applied to
# an active right-hand side, so an AD active type makes the spline value
# differentiable w.r.t. its knot values with no special handling.
#
# The XAD adjoint Tape<T,N> methods are instantiated only in src/Tape.cpp (in the
# odelia shared library), so this sourceCpp build links against it via PKG_LIBS,
# exactly like the leaf-thermal AD interface (see helper-load-odelia.R).

compile_spline_ad_interface <- function() {
  ensure_ode_interface_loaded(rebuild = FALSE)
  include_dir <- dirname(dirname(resolve_test_path(
    "include/odelia/ode_solver.hpp", "inst/include/odelia/ode_solver.hpp")))
  odelia_so <- .odelia_test_cache$odelia_so
  pkg_libs <- if (is.character(odelia_so) && length(odelia_so) == 1 &&
                  !is.na(odelia_so) && nzchar(odelia_so) && file.exists(odelia_so)) {
    shQuote(normalizePath(odelia_so, winslash = "/", mustWork = FALSE))
  } else {
    Sys.getenv("PKG_LIBS", unset = "")
  }
  withr::local_envvar(PKG_CPPFLAGS = paste0("-I", shQuote(include_dir)),
                      PKG_LIBS = pkg_libs)

  res <- tryCatch({
    Rcpp::sourceCpp(code = '
      #include <Rcpp.h>
      #include <vector>
      #include <XAD/XAD.hpp>
      #include <odelia/spline.hpp>

      // Reverse-mode gradient of the spline value at q w.r.t. each knot value.
      // [[Rcpp::export]]
      Rcpp::NumericVector spline_knot_value_grad(std::vector<double> x,
                                                 std::vector<double> y,
                                                 double q) {
        using ad = xad::adj<double>;
        using ad_t = ad::active_type;
        ad::tape_type tape;
        std::vector<ad_t> ya(y.begin(), y.end());
        for (auto& v : ya) tape.registerInput(v);
        tape.newRecording();
        odelia::spline::basic_spline<ad_t> s;
        s.set_points(x, ya);
        ad_t val = s(q);
        tape.registerOutput(val);
        xad::derivative(val) = 1.0;
        tape.computeAdjoints();
        Rcpp::NumericVector g(y.size());
        for (size_t i = 0; i < y.size(); ++i) g[i] = xad::derivative(ya[i]);
        return g;
      }

      // Double-path value through the same templated spline (alias path).
      // [[Rcpp::export]]
      double spline_value_double(std::vector<double> x, std::vector<double> y,
                                 double q) {
        odelia::spline::Spline s;   // = basic_spline<double>
        s.set_points(x, y);
        return s(q);
      }', verbose = FALSE)
    NULL
  }, error = function(e) e)

  if (inherits(res, "error")) {
    msg <- conditionMessage(res)
    if (grepl("active_tape_", msg, fixed = TRUE)) {
      testthat::skip("AD spline sourceCpp symbols unavailable in this load_all session; covered by installed-package tests.")
    }
    stop(res)
  }
}

testthat::test_that("basic_spline<ad> knot-value gradient matches finite differences", {
  testthat::skip_if(is_pkgload_dll(),
    "Skipping AD spline in pkgload load_all sessions due to native-pointer lifecycle.")
  compile_spline_ad_interface()

  x <- c(0.0, 1.0, 2.5, 4.0, 6.0, 9.0)
  y <- sin(0.4 * x) + 0.1 * x
  q <- 3.3

  g_ad <- spline_knot_value_grad(x, y, q)

  # Central FD of the value w.r.t. each knot value, through the double spline.
  h <- 1e-6
  g_fd <- vapply(seq_along(y), function(i) {
    yp <- y; ym <- y; yp[i] <- yp[i] + h; ym[i] <- ym[i] - h
    (spline_value_double(x, yp, q) - spline_value_double(x, ym, q)) / (2 * h)
  }, numeric(1))

  expect_equal(g_ad, g_fd, tolerance = 1e-6)

  # A cubic spline reproduces constants, so the sensitivities to all knot
  # values sum to exactly 1 at any query (partition of unity).
  expect_equal(sum(g_ad), 1, tolerance = 1e-9)
})

testthat::test_that("templated spline double alias is unchanged", {
  testthat::skip_if(is_pkgload_dll(), "Skipping sourceCpp build in load_all session.")
  compile_spline_ad_interface()

  # Cubic spline through samples of a cubic recovers it exactly.
  f <- function(x) 2 * x^3 - x^2 + 0.5 * x + 3
  x <- seq(-3, 3, length.out = 60)
  for (q in c(-2.1, 0.4, 1.7)) {
    expect_equal(spline_value_double(x, f(x), q), f(q), tolerance = 1e-8)
  }
})
