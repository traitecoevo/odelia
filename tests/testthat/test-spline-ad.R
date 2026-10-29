# Tests for the scalar-templated spline (odelia::interpolator::hermite_spline<S>), the
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
  # Without contraction, because read_linear is held to eval bit for bit: fused,
  # its end line's y + m * dx rounds once where the tape's arithmetic rounds twice.
  withr::local_envvar(PKG_CPPFLAGS = odelia_cppflags(include_dir),
                      PKG_CXXFLAGS = "-ffp-contract=off",
                      PKG_LIBS = pkg_libs)

  res <- tryCatch({
    Rcpp::sourceCpp(code = '
      #include <Rcpp.h>
      // [[Rcpp::plugins(cpp20)]]
      #include <vector>
      #include <XAD/XAD.hpp>
      #include <odelia/spline.hpp>
      #include <odelia/interpolator.hpp>
      #include <odelia/tangent.hpp>

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
        odelia::interpolator::hermite_spline<ad_t> s;
        s.init(x, ya, odelia::interpolator::natural_slopes<ad_t>(x, ya));
        ad_t val = s.eval(q);
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
        odelia::interpolator::hermite_spline<double> s;
        s.init(x, y, odelia::interpolator::natural_slopes<double>(x, y));
        return s.eval(q);
      }

      // Same, one level up through the Interpolator wrapper (what plant uses).
      // [[Rcpp::export]]
      Rcpp::NumericVector interp_knot_value_grad(std::vector<double> x,
                                                 std::vector<double> y,
                                                 double q) {
        using ad = xad::adj<double>;
        using ad_t = ad::active_type;
        ad::tape_type tape;
        std::vector<ad_t> ya(y.begin(), y.end());
        for (auto& v : ya) tape.registerInput(v);
        tape.newRecording();
        odelia::interpolator::hermite_interpolator<ad_t> I;
        I.init(x, ya);
        ad_t val = I.eval(q);
        tape.registerOutput(val);
        xad::derivative(val) = 1.0;
        tape.computeAdjoints();
        Rcpp::NumericVector g(y.size());
        for (size_t i = 0; i < y.size(); ++i) g[i] = xad::derivative(ya[i]);
        return g;
      }

      // A read at an ACTIVE position, with w as its adjoint: the value, its derivative
      // with respect to the position and to each knot value, and the statements
      // the read alone put on the tape. composed takes the value and slope at the
      // passive position and carries the query through the slope in a statement
      // of its own, which is the read eval makes outside the knots.
      // [[Rcpp::export]]
      Rcpp::List spline_active_query(std::vector<double> x, std::vector<double> y,
                                     double q, double w, bool composed) {
        using ad = xad::adj<double>;
        using ad_t = ad::active_type;
        ad::tape_type tape;
        std::vector<ad_t> ya(y.begin(), y.end());
        for (auto& v : ya) tape.registerInput(v);
        ad_t qa = q;
        tape.registerInput(qa);
        tape.newRecording();
        odelia::interpolator::hermite_spline<ad_t> s;
        s.init(x, ya, odelia::interpolator::natural_slopes<ad_t>(x, ya));
        const double before = tape.getNumStatements();
        ad_t val;
        if (composed) {
          ad_t v, dydu;
          s.value_and_slope(q, v, dydu);
          val = v + dydu * (qa - q);
        } else {
          val = s.eval(qa);
        }
        const double statements = tape.getNumStatements() - before;
        tape.registerOutput(val);
        xad::derivative(val) = w;
        tape.computeAdjoints();
        Rcpp::NumericVector g(y.size());
        for (size_t i = 0; i < y.size(); ++i) g[i] = xad::derivative(ya[i]);
        return Rcpp::List::create(
            Rcpp::_["value"] = xad::value(val),
            Rcpp::_["d_query"] = xad::derivative(qa),
            Rcpp::_["d_knots"] = g,
            Rcpp::_["statements"] = statements);
      }

      // [[Rcpp::export]]
      double interp_value_double(std::vector<double> x, std::vector<double> y,
                                 double q) {
        odelia::interpolator::Interpolator I;   // = hermite_interpolator<double>
        I.init(x, y);
        return I.eval(q);
      }

      // read_linear at q on a spline with active knot data: its value and slope
      // beside eval and slope, its weights summed per knot value and slope, and
      // the tape gradient of eval with respect to the same knot data.
      // [[Rcpp::export]]
      Rcpp::List spline_read_linear(std::vector<double> x, std::vector<double> y,
                                    std::vector<double> m, double q) {
        using ad = xad::adj<double>;
        using ad_t = ad::active_type;
        ad::tape_type tape;
        std::vector<ad_t> ya(y.begin(), y.end()), ma(m.begin(), m.end());
        for (auto& v : ya) tape.registerInput(v);
        for (auto& v : ma) tape.registerInput(v);
        tape.newRecording();
        odelia::interpolator::hermite_spline<ad_t> s;
        s.init(x, ya, ma);
        const auto r = s.read_linear(q);
        const ad_t* held_y = s.values().data();
        const ad_t* held_m = s.slopes().data();
        Rcpp::NumericVector wy(y.size()), wm(m.size());
        for (std::size_t i = 0; i < 4; ++i) {
          if (r.input[i] >= held_y && r.input[i] < held_y + y.size()) {
            wy[r.input[i] - held_y] += r.weight[i];
          } else {
            wm[r.input[i] - held_m] += r.weight[i];
          }
        }
        ad_t val = s.eval(q);
        tape.registerOutput(val);
        xad::derivative(val) = 1.0;
        tape.computeAdjoints();
        Rcpp::NumericVector gy(y.size()), gm(m.size());
        for (std::size_t i = 0; i < y.size(); ++i) {
          gy[i] = xad::derivative(ya[i]);
          gm[i] = xad::derivative(ma[i]);
        }
        return Rcpp::List::create(
            Rcpp::_["value"] = r.value, Rcpp::_["slope"] = r.slope,
            Rcpp::_["eval"] = xad::value(val),
            Rcpp::_["eval_slope"] = xad::value(s.slope(q)),
            Rcpp::_["weight_y"] = wy, Rcpp::_["weight_m"] = wm,
            Rcpp::_["grad_y"] = gy, Rcpp::_["grad_m"] = gm);
      }

      // a b / (a + b) and its two partials, from value_and_partials.
      // [[Rcpp::export]]
      Rcpp::NumericVector tangent_value_and_partials(double a, double b) {
        const auto r = odelia::ode::value_and_partials(
            [](const auto& p, const auto& q) { return p * q / (p + q); }, a, b);
        return Rcpp::NumericVector::create(r.value, r.partial[0], r.partial[1]);
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

testthat::test_that("hermite_spline<ad> knot-value gradient matches finite differences", {
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

testthat::test_that("hermite_interpolator<ad> matches FD (the wrapper plant uses)", {
  testthat::skip_if(is_pkgload_dll(),
    "Skipping AD interpolator in pkgload load_all sessions.")
  compile_spline_ad_interface()

  x <- c(0.0, 1.0, 2.5, 4.0, 6.0, 9.0)
  y <- sin(0.4 * x) + 0.1 * x
  q <- 3.3

  g_ad <- interp_knot_value_grad(x, y, q)
  h <- 1e-6
  g_fd <- vapply(seq_along(y), function(i) {
    yp <- y; ym <- y; yp[i] <- yp[i] + h; ym[i] <- ym[i] - h
    (interp_value_double(x, yp, q) - interp_value_double(x, ym, q)) / (2 * h)
  }, numeric(1))

  expect_equal(g_ad, g_fd, tolerance = 1e-6)
  expect_equal(sum(g_ad), 1, tolerance = 1e-9)
})

testthat::test_that("hermite_spline<ad> read at an active position", {
  testthat::skip_if(is_pkgload_dll(),
    "Skipping AD spline in pkgload load_all sessions due to native-pointer lifecycle.")
  compile_spline_ad_interface()

  x <- c(0.0, 1.0, 2.5, 4.0, 6.0, 9.0)
  y <- sin(0.4 * x) + 0.1 * x

  # Inside a span, on an interior knot, on each end knot, and past each end,
  # where the end span's line is extended. The adjoint is not one, because a unit
  # adjoint hides how a knot's partial is associated.
  w <- 0.3
  for (q in c(3.3, 0.37, 2.5, 0, 9, -0.5, 9.75)) {
    r <- spline_active_query(x, y, q, w, composed = FALSE)
    r0 <- spline_active_query(x, y, q, w, composed = TRUE)
    # The value and the position derivative are the composed read's bit for bit;
    # the knot adjoints agree to rounding.
    expect_identical(r$value, r0$value)
    expect_identical(r$d_query, r0$d_query)
    expect_equal(r$d_knots, r0$d_knots, tolerance = 1e-10)
    expect_equal(r$value, spline_value_double(x, y, q), tolerance = 1e-14)

    h <- 1e-6
    g_fd <- vapply(seq_along(y), function(i) {
      yp <- y; ym <- y; yp[i] <- yp[i] + h; ym[i] <- ym[i] - h
      (spline_value_double(x, yp, q) - spline_value_double(x, ym, q)) / (2 * h)
    }, numeric(1))
    expect_equal(r$d_knots, w * g_fd, tolerance = 1e-6)
  }

  # Inside the knots the read is one expression, so one statement, where the
  # composed read takes several.
  expect_equal(spline_active_query(x, y, 3.3, w, composed = FALSE)$statements, 1)
  expect_gt(spline_active_query(x, y, 3.3, w, composed = TRUE)$statements, 1)
})

testthat::test_that("read_linear is eval's read, and its weights are the read's gradient", {
  testthat::skip_if(is_pkgload_dll(), "Skipping sourceCpp build in load_all session.")
  compile_spline_ad_interface()

  x <- c(0.0, 1.0, 2.5, 4.0, 6.0, 9.0)
  y <- sin(0.4 * x) + 0.1 * x
  m <- 0.4 * cos(0.4 * x) + 0.1
  # Inside, on a knot, and beyond either end, where eval extends the end line.
  for (q in c(3.3, 0.2, 4.0, 8.99, -0.5, 9.7)) {
    r <- spline_read_linear(x, y, m, q)
    expect_identical(r$value, r$eval, label = paste("value at", q))
    expect_identical(r$slope, r$eval_slope, label = paste("slope at", q))
    # The weights are the Hermite basis; the tape reaches the knots through the
    # span coefficients, which associate differently. Absolute, because a weight
    # near zero (3e-5 at 8.99) has no relative precision to compare.
    expect_lt(max(abs(r$weight_y - r$grad_y)), 1e-14, label = paste("y weights at", q))
    expect_lt(max(abs(r$weight_m - r$grad_m)), 1e-14, label = paste("m weights at", q))
  }
})

testthat::test_that("value_and_partials is f and its exact partials", {
  testthat::skip_if(is_pkgload_dll(), "Skipping sourceCpp build in load_all session.")
  compile_spline_ad_interface()

  a <- 0.7; b <- 0.3
  r <- tangent_value_and_partials(a, b)
  expect_identical(r[1], a * b / (a + b))
  expect_equal(r[2:3], c(b^2, a^2) / (a + b)^2, tolerance = 1e-15)
})
