# Tests for odelia/interpolator.hpp: the C1 piecewise cubic carrying a value and a
# slope at each knot, the slope rule for knots that arrive with values alone, and
# the node rule for a target whose features are not known in advance.
#
# The interpolant records on a tape when its scalar or its query position is active,
# and the XAD Tape<T,N> template methods are explicitly instantiated only in the
# odelia shared library, so these snippets link against it rather than compiling
# header-only.

compile_interpolator_interface <- function() {
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
    #include <odelia/interpolator.hpp>

    using tape_type = xad::Tape<double>;
    using adouble = tape_type::active_type;
    using odelia::interpolator::hermite_spline;

    // A vertical field sampled at knots, with the slope from the same expression the
    // value comes from -- the pair a crown integral reads.
    struct AttenuatingField {
      double extinction = 1.7, amplitude = 2.3;

      double value(double z) const { return amplitude * std::exp(-extinction * z); }
      double slope(double z) const { return -extinction * value(z); }
      double curvature(double z) const { return extinction * extinction * value(z); }
    };

    static AttenuatingField field;

    static std::vector<double> field_values(const std::vector<double>& z) {
      std::vector<double> v(z.size());
      for (std::size_t i = 0; i < z.size(); ++i) v[i] = field.value(z[i]);
      return v;
    }
    static std::vector<double> field_slopes(const std::vector<double>& z) {
      std::vector<double> v(z.size());
      for (std::size_t i = 0; i < z.size(); ++i) v[i] = field.slope(z[i]);
      return v;
    }
    static std::vector<double> field_curvatures(const std::vector<double>& z) {
      std::vector<double> v(z.size());
      for (std::size_t i = 0; i < z.size(); ++i) v[i] = field.curvature(z[i]);
      return v;
    }

    // [[Rcpp::export]]
    Rcpp::List hermite_value_and_slope(std::vector<double> z, std::vector<double> u) {
      hermite_spline<double> interp;
      interp.set_nodes(z);
      interp.set_data(field_values(z), field_slopes(z));
      std::vector<double> value(u.size()), slope(u.size()),
                          pair_value(u.size()), pair_slope(u.size());
      for (std::size_t i = 0; i < u.size(); ++i) {
        value[i] = interp.eval(u[i]);
        slope[i] = interp.slope(u[i]);
        interp.value_and_slope(u[i], pair_value[i], pair_slope[i]);
      }
      return Rcpp::List::create(Rcpp::_["value"] = value, Rcpp::_["slope"] = slope,
                                Rcpp::_["pair_value"] = pair_value,
                                Rcpp::_["pair_slope"] = pair_slope);
    }

    // Setting the nodes once and the data twice, against two all-at-once builds on the
    // same nodes: what a per-stage rebuild does, compared with rebuilding everything.
    // [[Rcpp::export]]
    Rcpp::List hermite_split_matches_all_at_once(std::vector<double> z,
                                                 std::vector<double> u,
                                                 double scale) {
      const std::vector<double> y = field_values(z), m = field_slopes(z);
      std::vector<double> y2(y.size()), m2(m.size());
      for (std::size_t i = 0; i < y.size(); ++i) { y2[i] = scale * y[i]; m2[i] = scale * m[i]; }

      hermite_spline<double> split, whole;
      split.set_nodes(z);

      std::vector<int> same(2 * u.size());
      split.set_data(y, m);
      whole.init(z, y, m);
      for (std::size_t i = 0; i < u.size(); ++i) {
        same[i] = odelia::util::identical(split.eval(u[i]), whole.eval(u[i])) &&
                  odelia::util::identical(split.slope(u[i]), whole.slope(u[i]));
      }
      // The second set_data reuses the span layout the first one wrote into.
      split.set_data(y2, m2);
      whole.init(z, y2, m2);
      for (std::size_t i = 0; i < u.size(); ++i) {
        same[u.size() + i] =
            odelia::util::identical(split.eval(u[i]), whole.eval(u[i])) &&
            odelia::util::identical(split.slope(u[i]), whole.slope(u[i]));
      }
      return Rcpp::List::create(Rcpp::_["identical"] = same,
                                Rcpp::_["initialised"] = split.is_initialised(),
                                Rcpp::_["size"] = static_cast<int>(split.size()));
    }

    // d(value)/d(u) recorded by a read at an active position.
    // [[Rcpp::export]]
    std::vector<double> hermite_active_position_adjoint(std::vector<double> z,
                                                        std::vector<double> u) {
      const std::vector<double> yd = field_values(z), md = field_slopes(z);
      const std::vector<adouble> y(yd.begin(), yd.end()), m(md.begin(), md.end());
      hermite_spline<adouble> interp;
      interp.set_nodes(z);
      interp.set_data(y, m);

      std::vector<double> d(u.size());
      for (std::size_t i = 0; i < u.size(); ++i) {
        tape_type tape;
        adouble ui = u[i];
        tape.registerInput(ui);
        tape.newRecording();
        adouble value = interp.eval(ui);
        tape.registerOutput(value);
        xad::derivative(value) = 1.0;
        tape.computeAdjoints();
        d[i] = xad::derivative(ui);
      }
      return d;
    }

    // [[Rcpp::export]]
    double hermite_eval_before_data(std::vector<double> z) {
      hermite_spline<double> interp;
      interp.set_nodes(z);
      return interp.eval(z.front());
    }

    // [[Rcpp::export]]
    void hermite_data_before_nodes(std::vector<double> z) {
      hermite_spline<double> interp;
      interp.set_data(field_values(z), field_slopes(z));
    }

    // [[Rcpp::export]]
    void hermite_set_nodes_descending(std::vector<double> z) {
      hermite_spline<double> interp;
      interp.set_nodes(z);
    }

    // A cubic given its own values and slopes is reproduced exactly, which is what
    // says the coefficients belong to that polynomial and not to a fit.
    // [[Rcpp::export]]
    std::vector<double> hermite_cubic_exact(std::vector<double> z,
                                            std::vector<double> u) {
      auto f  = [](double t) { return 0.7 + 1.3 * t - 0.4 * t * t + 0.11 * t * t * t; };
      auto df = [](double t) { return 1.3 - 0.8 * t + 0.33 * t * t; };
      std::vector<double> y(z.size()), m(z.size());
      for (std::size_t i = 0; i < z.size(); ++i) { y[i] = f(z[i]); m[i] = df(z[i]); }
      hermite_spline<double> interp;
      interp.init(z, y, m);
      std::vector<double> out(u.size());
      for (std::size_t i = 0; i < u.size(); ++i) out[i] = interp.eval(u[i]) - f(u[i]);
      return out;
    }

    // The adjoint of one read with respect to every knot VALUE. Each span reads two
    // knots, so a read touches four of the 2K inputs and the rest are exactly zero.
    // [[Rcpp::export]]
    std::vector<double> hermite_knot_value_adjoint(std::vector<double> z, double u) {
      tape_type tape;
      const std::vector<double> y0 = field_values(z), m0 = field_slopes(z);
      std::vector<adouble> y(y0.begin(), y0.end()), m(m0.begin(), m0.end());
      tape.registerInputs(y);
      tape.registerInputs(m);
      tape.newRecording();
      hermite_spline<adouble> interp;
      interp.init(z, y, m);
      adouble out = interp.eval(u);
      tape.registerOutput(out);
      xad::derivative(out) = 1.0;
      tape.computeAdjoints();
      std::vector<double> d;
      d.reserve(y.size() + m.size());
      for (const adouble& v : y) d.push_back(xad::derivative(v));
      for (const adouble& v : m) d.push_back(xad::derivative(v));
      return d;
    }

    // [[Rcpp::export]]
    std::vector<double> interpolator_monotone_slopes(std::vector<double> x,
                                                     std::vector<double> y) {
      return odelia::interpolator::monotone_slopes(x, y);
    }

    // Read a limited fit between its control points, to see whether it leaves their
    // range.
    // [[Rcpp::export]]
    std::vector<double> interpolator_monotone_read(std::vector<double> x,
                                                   std::vector<double> y,
                                                   std::vector<double> u) {
      hermite_spline<double> interp;
      interp.init(x, y, odelia::interpolator::monotone_slopes(x, y));
      std::vector<double> out(u.size());
      for (std::size_t i = 0; i < u.size(); ++i) out[i] = interp.eval(u[i]);
      return out;
    }

    // The nodes refinement places on the field, and how well the fit on them holds.
    // [[Rcpp::export]]
    Rcpp::List interpolator_refine(double a, double b, double tol) {
      const auto chosen = odelia::interpolator::refine<double>(
          [](double z) { return std::pair<double, double>(field.value(z),
                                                          field.slope(z)); },
          a, b, tol);
      hermite_spline<double> interp;
      interp.init(chosen.x, chosen.y, chosen.m);
      double worst = 0.0;
      for (std::size_t k = 0; k + 1 < chosen.x.size(); ++k) {
        const double mid = 0.5 * (chosen.x[k] + chosen.x[k + 1]);
        worst = std::max(worst, std::abs(interp.eval(mid) - field.value(mid)));
      }
      return Rcpp::List::create(Rcpp::_["n"] = (int)chosen.x.size(),
                                Rcpp::_["worst"] = worst,
                                Rcpp::_["x"] = chosen.x);
    }

    // The same knots read at both orders, and each order at its own knots. A source
    // with a closed form for the curvature can supply three channels; one with only
    // two cannot, and set_data has a signature for each.
    // [[Rcpp::export]]
    Rcpp::List hermite_orders(std::vector<double> z, std::vector<double> u) {
      hermite_spline<double> cubic;
      cubic.init(z, field_values(z), field_slopes(z));
      hermite_spline<double, 5> quintic;
      quintic.init(z, field_values(z), field_slopes(z), field_curvatures(z));
      double worst3 = 0.0, worst5 = 0.0, at_knots = 0.0, slope_at_knots = 0.0;
      for (std::size_t i = 0; i < u.size(); ++i) {
        worst3 = std::max(worst3, std::abs(cubic.eval(u[i]) - field.value(u[i])));
        worst5 = std::max(worst5, std::abs(quintic.eval(u[i]) - field.value(u[i])));
      }
      for (std::size_t k = 0; k < z.size(); ++k) {
        at_knots = std::max(at_knots, std::abs(quintic.eval(z[k]) - field.value(z[k])));
        slope_at_knots =
            std::max(slope_at_knots, std::abs(quintic.slope(z[k]) - field.slope(z[k])));
      }
      // The slope against a difference of the quintic\'s OWN value, so this checks
      // one polynomial rather than agreement with the field.
      const double h = 1e-6;
      double own = 0.0;
      for (std::size_t i = 0; i < u.size(); ++i) {
        if (u[i] - h <= z.front() || u[i] + h >= z.back()) continue;
        own = std::max(own, std::abs((quintic.eval(u[i] + h) - quintic.eval(u[i] - h))
                                     / (2 * h) - quintic.slope(u[i])));
      }
      return Rcpp::List::create(Rcpp::_["cubic"] = worst3,
                                Rcpp::_["quintic"] = worst5,
                                Rcpp::_["at_knots"] = at_knots,
                                Rcpp::_["slope_at_knots"] = slope_at_knots,
                                Rcpp::_["own_slope"] = own);
    }', verbose = FALSE)
}

# Two knot sets: uniform (indexed by arithmetic) and not (indexed by a search).
knot_sets <- list(uniform = seq(0, 4, length.out = 33),
                  irregular = c(0, 0.07, 0.31, 0.4, 0.95, 1.6, 2.05, 2.9, 3.4, 4))
query <- c(0.05, 0.31, 0.62, 1.0, 1.55, 2.4, 3.87, 3.999)

testthat::test_that("slope is the exact derivative of what eval returns", {
  compile_interpolator_interface()

  # The slope is not a fit by-product: it is the derivative of the same cubic eval
  # uses, so it agrees with a central difference of eval everywhere, knots included.
  h <- 1e-6
  for (z in knot_sets) {
    r <- hermite_value_and_slope(z, query)
    fd <- (hermite_value_and_slope(z, query + h)$value -
           hermite_value_and_slope(z, query - h)$value) / (2 * h)
    testthat::expect_equal(r$slope, fd, tolerance = 1e-8)

    # value_and_slope loads the span once; it must return the same pair.
    testthat::expect_identical(r$pair_value, r$value)
    testthat::expect_identical(r$pair_slope, r$slope)
  }
})

testthat::test_that("setting the nodes once and the data per stage moves no value", {
  compile_interpolator_interface()

  for (z in knot_sets) {
    r <- hermite_split_matches_all_at_once(z, query, scale = 0.37)
    testthat::expect_true(all(r$identical == 1L))
    testthat::expect_true(r$initialised)
    testthat::expect_equal(r$size, length(z))
  }
})

testthat::test_that("a read at an active position records the derivative of the query", {
  compile_interpolator_interface()

  # Without the with_query_derivative the span is indexed at the passive position and nothing records
  # the query, so this adjoint would be exactly zero.
  h <- 1e-6
  for (z in knot_sets) {
    adjoint <- hermite_active_position_adjoint(z, query)
    testthat::expect_false(all(adjoint == 0))

    fd <- (hermite_value_and_slope(z, query + h)$value -
           hermite_value_and_slope(z, query - h)$value) / (2 * h)
    worst <- max(abs(adjoint - fd) / abs(fd))
    message(sprintf("active-position adjoint vs central difference (h=%g): worst rel=%.3e",
                    h, worst))
    testthat::expect_lt(worst, 1e-7)
  }
})

testthat::test_that("an incomplete or unusable build stops", {
  compile_interpolator_interface()

  z <- knot_sets$uniform
  testthat::expect_error(hermite_eval_before_data(z), "not initialised")
  testthat::expect_error(hermite_data_before_nodes(z), "no knots set")
  testthat::expect_error(hermite_set_nodes_descending(rev(z)), "strictly ascending")
  testthat::expect_error(hermite_set_nodes_descending(c(1.0)), "at least 2 knots")
})

testthat::test_that("a cubic given its own values and slopes comes back exactly", {
  compile_interpolator_interface()
  for (z in knot_sets) {
    testthat::expect_lt(max(abs(hermite_cubic_exact(z, query))), 1e-13)
  }
})

testthat::test_that("a read touches four knot inputs and no others", {
  compile_interpolator_interface()
  z <- knot_sets$uniform
  # A query strictly inside one span: the two knots bounding it, each through its
  # value and its slope. Everything else must be exactly zero, not merely small --
  # this is the sparsity an O(1) adjoint rests on.
  d <- hermite_knot_value_adjoint(z, 1.55)
  testthat::expect_equal(sum(d != 0), 4L)
  k <- findInterval(1.55, z)
  nz <- c(k, k + 1L, length(z) + k, length(z) + k + 1L)
  testthat::expect_equal(which(d != 0), nz)
})

testthat::test_that("the slope rule keeps a limited fit inside its own data", {
  compile_interpolator_interface()
  # An intermittent series: mostly zero, occasional pulses. A fit that chooses
  # slopes globally reads negative between two dry points; this one cannot.
  x <- seq(0, 20, length.out = 41)
  y <- ifelse(seq_along(x) %% 7 == 0, 5, 0)
  u <- seq(0, 20, length.out = 1000)
  got <- interpolator_monotone_read(x, y, u)
  testthat::expect_gte(min(got), 0)
  testthat::expect_lte(max(got), max(y))

  # And on a monotone series it is monotone.
  ym <- cumsum(c(0, abs(sin(seq_along(x)[-1]))))
  gotm <- interpolator_monotone_read(x, ym, u)
  testthat::expect_true(all(diff(gotm) >= -1e-12))

  # A flat pair pins both its slopes to zero, which is what stops a pulse being
  # smeared backwards into a dry day.
  m <- interpolator_monotone_slopes(x, y)
  testthat::expect_equal(m[2], 0)
})

testthat::test_that("refinement stops where the fit resolves the target", {
  compile_interpolator_interface()
  loose <- interpolator_refine(0, 4, 1e-4)
  tight <- interpolator_refine(0, 4, 1e-8)
  # A tighter tolerance places more nodes and reaches a smaller miss.
  testthat::expect_gt(tight$n, loose$n)
  testthat::expect_lt(tight$worst, loose$worst)
  testthat::expect_lt(loose$worst, 1e-4)
  # The nodes it placed are a usable knot set: ascending, spanning the interval.
  testthat::expect_true(all(diff(loose$x) > 0))
  testthat::expect_equal(range(loose$x), c(0, 4))
})

testthat::test_that("a curvature channel raises the order and is refused where absent", {
  compile_interpolator_interface()

  # Both orders on the same knots, read through the one polynomial the class writes
  # for both. A cubic is determined by a value and a slope, a quintic by those plus a
  # curvature, so the second converges as h^6 against h^4 -- which is what lets a
  # source with a closed-form second derivative reach a given error on far fewer
  # knots. Measured here on a coarse grid so the two are apart by more than the
  # arithmetic's own noise.
  coarse <- seq(0, 4, length.out = 9)
  probe <- seq(0.02, 3.98, length.out = 401)
  got <- hermite_orders(coarse, probe)

  testthat::expect_lt(got$quintic, got$cubic / 100)

  # Both channels are exact at the knots the data was supplied at: a quintic is an
  # interpolant, not a fit, in the value and the slope alike.
  testthat::expect_identical(got$at_knots, 0)
  testthat::expect_identical(got$slope_at_knots, 0)

  # And its slope is the derivative of its own value, not of the target -- the same
  # property the cubic has, checked against the polynomial rather than the field.
  testthat::expect_lt(got$own_slope, 1e-8)
})
