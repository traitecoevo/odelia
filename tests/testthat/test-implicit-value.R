# Tests for the two ways a value computed away from the tape gets onto it:
# record_with_derivatives, which takes the derivatives as supplied numbers, and
# implicit_value, which derives them from a residual through the implicit
# function theorem and records the result through the first.
#
# Both record on a tape, and the XAD Tape<T,N> template methods are
# explicitly instantiated only in the odelia shared library, so these snippets link
# against it rather than compiling header-only.

compile_implicit_value_interface <- function() {
  ensure_ode_interface_loaded()

  include_dir <- odelia_include_dir()
  odelia_so <- .odelia_test_cache$odelia_so
  withr::local_envvar(
    PKG_CPPFLAGS = odelia_cppflags(include_dir),
    PKG_LIBS = shQuote(normalizePath(odelia_so, winslash = "/", mustWork = TRUE))
  )
  Rcpp::sourceCpp(code = '
    // [[Rcpp::plugins(cpp20)]]
    #include <Rcpp.h>
    #include <vector>
    #include <cmath>
    #include <XAD/XAD.hpp>
    #include <odelia/implicit_node.hpp>

    using tape_type = xad::Tape<double>;
    using adouble = tape_type::active_type;

    // A two-parameter system whose operating point is defined implicitly: the root of
    // a*y^3 + y - b = 0. Written once in the working scalar, so the double root-find
    // and the recorded evaluation read the same equation.
    template <typename T>
    struct CubicBalance {
      using value_type = T;
      T a, b;

      // Declared -> T: a deduced return type here is an XAD expression template
      // holding references to the temporaries of this return statement.
      T residual(T y) const { return a * y * y * y + y - b; }

      // The operating point, in double and off the tape.
      double solve() const {
        double y = 0.0;
        for (int i = 0; i < 100; ++i) {
          const double ad = odelia::util::to_passive(a);
          const double bd = odelia::util::to_passive(b);
          const double f = ad * y * y * y + y - bd;
          const double fp = 3.0 * ad * y * y + 1.0;
          const double dy = f / fp;
          y -= dy;
          if (std::abs(dy) < 1e-15 * (std::abs(y) + 1.0)) break;
        }
        return y;
      }
    };

    // [[Rcpp::export]]
    double cubic_root(double a, double b) {
      CubicBalance<double> sys{a, b};
      return sys.solve();
    }

    // The value and d(y*)/d(a), d(y*)/d(b) from one reverse sweep.
    // [[Rcpp::export]]
    Rcpp::List implicit_value_gradient(double a, double b) {
      tape_type tape;
      CubicBalance<adouble> sys{adouble(a), adouble(b)};
      tape.registerInput(sys.a);
      tape.registerInput(sys.b);
      tape.newRecording();

      const double y_star = sys.solve();
      // d/dy (a y^3 + y - b) = 3 a y^2 + 1, at the point the theorem divides at.
      adouble y = odelia::implicit_value<adouble>(
          y_star, 3.0 * a * y_star * y_star + 1.0,
          [&](const adouble& yy) -> adouble { return sys.residual(yy); });

      tape.registerOutput(y);
      xad::derivative(y) = 1.0;
      tape.computeAdjoints();
      return Rcpp::List::create(
          Rcpp::_["value"] = xad::value(y),
          Rcpp::_["y_star"] = y_star,
          Rcpp::_["identical"] = odelia::util::identical(xad::value(y), y_star),
          Rcpp::_["d_da"] = xad::derivative(sys.a),
          Rcpp::_["d_db"] = xad::derivative(sys.b));
    }

    // A value the tape never saw computed, entering it against two inputs it
    // never touched. Nothing here relates `value` to u or v, so the derivatives
    // that come back can only be the supplied ones.
    // [[Rcpp::export]]
    Rcpp::List record_with_derivatives_gradient(double value, double du,
                                                double dv) {
      tape_type tape;
      adouble u(1.5), v(-0.25);
      tape.registerInput(u);
      tape.registerInput(v);
      tape.newRecording();

      adouble y;
      // A named container, not a braced list: std::span has no initializer_list
      // constructor before C++26, so `{{u, du}, {v, dv}}` compiles only on
      // libstdc++ (which ships that constructor early) and not on libc++.
      const std::vector<odelia::input_and_derivative<adouble>> rows{{u, du},
                                                                    {v, dv}};
      const odelia::record_report report =
          odelia::record_with_derivatives<adouble>(value, rows, y);
      tape.registerOutput(y);
      xad::derivative(y) = 1.0;
      tape.computeAdjoints();
      return Rcpp::List::create(
          Rcpp::_["value"] = xad::value(y),
          Rcpp::_["identical"] = odelia::util::identical(xad::value(y), value),
          Rcpp::_["d_du"] = xad::derivative(u),
          Rcpp::_["d_dv"] = xad::derivative(v),
          Rcpp::_["whole"] = report.whole,
          Rcpp::_["at"] = (int) report.at,
          Rcpp::_["why"] = report.why);
    }

    // A root p of R(p; u, v) = 0 and one output that depends on p recorded against
    // it. The residual is written so that its derivatives ARE the two slopes and its
    // value at the root is exactly zero, so what comes back can only be the quotient
    // the theorem gives and the chain through p.
    //
    // Every kind closes on its residual, so this exercises the theorem through
    // implicit_value rather than through supplied slope rows.
    // [[Rcpp::export]]
    Rcpp::List implicit_root_gradient(double p, double residual_slope,
                                      double dR_du, double dR_dv, double dy_dp) {
      tape_type tape;
      adouble u(1.5), v(-0.25);
      tape.registerInput(u);
      tape.registerInput(v);
      tape.newRecording();

      const double u0 = 1.5, v0 = -0.25;
      adouble root;
      odelia::record_report on_root{true, 0, ""};
      try {
        root = odelia::implicit_value<adouble>(
            p, residual_slope, [&](const adouble& y) -> adouble {
              (void)y;
              return adouble(dR_du * (u - u0) + dR_dv * (v - v0));
            });
      } catch (const std::runtime_error& e) {
        // A fold and a non-finite row both stop here now: the reporting sibling
        // went because its one consumer turned the report into the same refusal
        // this exception becomes. Caught so the checks below can read which.
        root = adouble(p);
        on_root = odelia::record_report{false, 0, e.what()};
      }
      adouble y;
      const std::vector<odelia::input_and_derivative<adouble>> on_root_rows{
          {root, dy_dp}};
      const odelia::record_report on_y =
          odelia::record_with_derivatives<adouble>(7.5, on_root_rows, y);
      tape.registerOutput(y);
      xad::derivative(y) = 1.0;
      tape.computeAdjoints();
      return Rcpp::List::create(
          Rcpp::_["root"] = xad::value(root),
          Rcpp::_["identical"] = odelia::util::identical(xad::value(root), p),
          Rcpp::_["y"] = xad::value(y),
          Rcpp::_["y_identical"] = odelia::util::identical(xad::value(y), 7.5),
          Rcpp::_["d_du"] = xad::derivative(u),
          Rcpp::_["d_dv"] = xad::derivative(v),
          Rcpp::_["root_whole"] = on_root.whole,
          Rcpp::_["root_why"] = on_root.why,
          Rcpp::_["y_whole"] = on_y.whole,
          Rcpp::_["y_why"] = on_y.why);
    }

    // A residual that touches zero rather than crossing it, so dF/dy is zero at
    // the operating point and the quotient the theorem asks for does not exist.
    // [[Rcpp::export]]
    double implicit_value_at_fold(double a) {
      tape_type tape;
      adouble aa(a);
      tape.registerInput(aa);
      tape.newRecording();
      // d/dy (y - a)^2 = 2(y - a), which is exactly zero at y = a: the fold.
      adouble y = odelia::implicit_value<adouble>(
          a, 0.0,
          [&](const adouble& yy) -> adouble { return (yy - aa) * (yy - aa); });
      return xad::value(y);
    }', verbose = FALSE)
}

testthat::test_that("record_with_derivatives returns the value and carries what it was given", {
  compile_implicit_value_interface()

  got <- record_with_derivatives_gradient(4.75, 2.0, -3.5)

  # Each term is an input minus its own passive copy, which is zero in value, so
  # what comes back is the number handed in rather than one near it.
  testthat::expect_true(got$identical)
  testthat::expect_equal(got$d_du, 2.0)
  testthat::expect_equal(got$d_dv, -3.5)

  testthat::expect_true(got$whole)

  # A non-finite derivative poisons the VALUE and not only the tape, because NaN
  # times zero is not a number. So the rows are refused -- ALL of them, not the
  # one that was bad: a value carrying part of its rows is a channel gone missing
  # with every number still finite. The VALUE is handed back either way, because
  # whether a caller can go on without the rows is the caller's to decide.
  bad <- record_with_derivatives_gradient(4.75, NaN, 1.0)
  testthat::expect_false(bad$whole)
  testthat::expect_match(bad$why, "is not finite")
  testthat::expect_true(bad$identical)
  testthat::expect_equal(bad$at, 0L)
  testthat::expect_equal(bad$d_du, 0.0)
  testthat::expect_equal(bad$d_dv, 0.0)
})

testthat::test_that("implicit_value returns the operating point and its IFT derivative", {
  compile_implicit_value_interface()

  a <- 0.35
  b <- 2.2
  got <- implicit_value_gradient(a, b)

  # The node returns y* itself, so a parameter the equation does not reach can
  # introduce no shift.
  testthat::expect_true(got$identical)

  # Central difference of the same root-find, at three step sizes. The node's dF/dy
  # is itself a central difference at eps = 1e-6*(|y*|+1), so this reports what that
  # costs rather than asserting it away.
  eps <- c(1e-3, 1e-4, 1e-5)
  fd <- vapply(eps, function(h) {
    c((cubic_root(a + h, b) - cubic_root(a - h, b)) / (2 * h),
      (cubic_root(a, b + h) - cubic_root(a, b - h)) / (2 * h))
  }, numeric(2))

  residual <- vapply(seq_along(eps), function(i) {
    max(abs(c(got$d_da, got$d_db) - fd[, i]) / abs(fd[, i]))
  }, numeric(1))
  message(sprintf("IFT vs central difference: eps=%g rel=%.3e | eps=%g rel=%.3e | eps=%g rel=%.3e",
                  eps[[1]], residual[[1]], eps[[2]], residual[[2]],
                  eps[[3]], residual[[3]]))

  # The residual falls as eps^2 across the three, so what it measures is the
  # reference difference's own truncation, not the node's dF/dy: the node's error is
  # below the smallest of these.
  for (i in seq_along(eps)) {
    testthat::expect_lt(residual[[i]], 5 * eps[[i]]^2)
  }
})

testthat::test_that("the theorem turns a residual's rows into its quotient's quotient", {
  compile_implicit_value_interface()

  p <- -1.75
  slope <- -0.8
  got <- implicit_root_gradient(p, slope, 2.0, -3.5, 1.0)

  # The root comes back as the number the solve left, so an input the residual
  # does not reach can introduce no shift.
  testthat::expect_true(got$identical)
  testthat::expect_equal(got$d_du, -2.0 / slope)
  testthat::expect_equal(got$d_dv, 3.5 / slope)

  # An output recorded against the root picks the quotient up through its own
  # slope, which is the whole reason the root is a value rather than a table of
  # per-input quotients: one root, any number of outputs.
  chained <- implicit_root_gradient(p, slope, 2.0, -3.5, 4.25)
  testthat::expect_true(chained$y_identical)
  testthat::expect_equal(chained$d_du, 4.25 * (-2.0 / slope))
  testthat::expect_equal(chained$d_dv, 4.25 * (3.5 / slope))
})

testthat::test_that("the reporting form says where the theorem does not apply", {
  compile_implicit_value_interface()

  # dR/dp of zero is the fold: the quotient is garbage rather than large, and a
  # non-finite one has nothing to divide by at all. Reported rather than thrown,
  # because the point is still the point and a consumer's other outputs may not
  # read it -- a stop takes those too.
  for (slope in c(0.0, NaN)) {
    at_fold <- implicit_root_gradient(1.0, slope, 2.0, -3.5, 1.0)
    testthat::expect_false(at_fold$root_whole)
    testthat::expect_match(at_fold$root_why, "does not apply")
    testthat::expect_true(at_fold$identical)
    testthat::expect_equal(at_fold$d_du, 0.0)
  }

  # And a residual row that is not finite is still refused by the record it records
  # through, after the quotient rather than before it.
  bad_slope <- implicit_root_gradient(1.0, -0.8, NaN, -3.5, 1.0)
  testthat::expect_false(bad_slope$root_whole)
  testthat::expect_match(bad_slope$root_why, "is not finite")
})

testthat::test_that("implicit_value stops where the theorem does not apply", {
  compile_implicit_value_interface()

  # F(y) = (y - a)^2 is zero at y* = a and so is its slope there, so the node must
  # stop rather than divide by it and return a gradient nothing supports.
  testthat::expect_error(implicit_value_at_fold(0.35), "does not apply")
})

