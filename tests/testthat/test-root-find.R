# root_find.hpp carries Boost.Math's TOMS748 and bisection so that the header
# core needs no Boost. These tests run both side by side, on the test problems
# of Alefeld, Potra and Shi (1995) at several tolerances and iteration caps, and
# require the same search: every point the target is evaluated at, the final
# bracket and the iterations spent, bit for bit. A cap small enough that the
# search runs out is included, since that path is where the two differ in what
# they report: Boost returns the bracket, odelia refuses it.
#
# The contract itself (refusals, exact zeros, no allocation) is tested in plain
# C++ in tests/standalone/root_find.cpp.

compile_root_find_comparison <- function() {
  if (isTRUE(.odelia_test_cache$root_find_loaded)) {
    return(invisible())
  }
  withr::local_envvar(PKG_CPPFLAGS = odelia_cppflags())
  Rcpp::sourceCpp(code = '
    // [[Rcpp::plugins(cpp20)]]
    // [[Rcpp::depends(BH)]]
    #include <Rcpp.h>
    #include <cmath>
    #include <vector>
    #include <boost/math/tools/roots.hpp>
    #include <boost/math/tools/toms748_solve.hpp>
    #include <odelia/root_find.hpp>

    // Alefeld, Potra and Shi (1995), Table I: problem `k` at parameter `n`.
    double alefeld(int k, double n, double x) {
      switch (k) {
      case 1: return std::sin(x) - x / 2;
      case 2: {
        double s = 0;
        for (int i = 1; i <= 20; ++i) {
          s += (2 * i - 5) * (2 * i - 5) / std::pow(x - i * i, 3);
        }
        return -2 * s;
      }
      case 3: return n * x * std::exp(-x * std::abs(n) / 40);
      case 4: return std::pow(x, n) - 0.2;
      case 5: return std::sin(x) - 0.5;
      case 6: return 2 * x * std::exp(-n) - 2 * std::exp(-n * x) + 1;
      case 7: return (1 + (1 - n) * (1 - n)) * x - (1 - n * x) * (1 - n * x);
      case 8: return x * x - std::pow(1 - x, n);
      case 9: return (1 + std::pow(1 - n, 4)) * x - std::pow(1 - n * x, 4);
      case 10: return std::exp(-n * x) * (x - 1) + std::pow(x, n);
      case 11: return (n * x - 1) / ((n - 1) * x);
      case 12: return std::pow(x, 1 / n) - std::pow(n, 1 / n);
      case 13: return x == 0 ? 0 : x / std::exp(1 / (x * x));
      case 14: return x >= 0 ? (n / 20) * (x / 1.5 + std::sin(x) - 1) : -n / 20;
      case 15: {
        const double e = 1.859;
        if (x > 2e-3 / (1 + n)) return std::exp(1.0) - e;
        if (x < 0) return -0.859;
        return std::exp((n + 1) * x / 2 * 1000) - e;
      }
      default: return x < 0.3 ? -1.0 : 1.0;  // a step
      }
    }

    struct recorded {
      int k;
      double n;
      std::vector<double>* xs;
      double operator()(double x) const {
        xs->push_back(x);
        return alefeld(k, n, x);
      }
    };

    // One search by each of the two, as lists of what each did. `form` is
    // "bracket" (TOMS748 handed the ends\' values), "ends" (TOMS748 evaluating
    // them) or "bisect". A Boost search that ran out is reported, not refused,
    // so whether it met the tolerance is returned for the test to read.
    // [[Rcpp::export]]
    Rcpp::List root_find_both(std::string form, int k, double n, double a,
                              double b, double atol, double rtol, int cap) {
      std::vector<double> boost_xs, odelia_xs;
      recorded fb{k, n, &boost_xs}, fo{k, n, &odelia_xs};
      odelia::root::tolerance tol{atol, rtol, static_cast<std::size_t>(cap)};
      auto boost_tol = [&tol](double x, double y) { return tol.met(x, y); };
      boost::uintmax_t boost_it = cap;
      std::uintmax_t odelia_it = cap;
      std::pair<double, double> rb(NA_REAL, NA_REAL), ro(NA_REAL, NA_REAL);
      std::string boost_err = "", odelia_err = "";
      double odelia_mid = NA_REAL;
      try {
        if (form == "bracket") {
          rb = boost::math::tools::toms748_solve(fb, a, b, alefeld(k, n, a),
                                                 alefeld(k, n, b), boost_tol,
                                                 boost_it);
        } else if (form == "ends") {
          rb = boost::math::tools::toms748_solve(fb, a, b, boost_tol, boost_it);
        } else {
          rb = boost::math::tools::bisect(fb, a, b, boost_tol, boost_it);
        }
      } catch (const std::exception& e) {
        boost_err = e.what();
      }
      try {
        if (form == "bracket") {
          ro = odelia::root::internal::toms748_solve(
              fo, a, b, alefeld(k, n, a), alefeld(k, n, b), tol, odelia_it);
          const auto rp = odelia::root::toms748_bracket(
              [&](double x) { return alefeld(k, n, x); }, a, b,
              alefeld(k, n, a), alefeld(k, n, b), tol);
          if (rp != ro) {
            odelia_err = "toms748_bracket differs from the search";
          }
          odelia_mid = odelia::root::toms748(
              [&](double x) { return alefeld(k, n, x); }, a, b,
              alefeld(k, n, a), alefeld(k, n, b), tol);
        } else if (form == "ends") {
          odelia_mid = odelia::root::toms748(fo, a, b, tol);
        } else {
          odelia_mid = odelia::root::bisect(fo, a, b, tol);
        }
      } catch (const odelia::util::RootFindError& e) {
        odelia_err = std::string("RootFindError: ") + e.what();
      } catch (const odelia::util::BracketError& e) {
        odelia_err = std::string("BracketError: ") + e.what();
      }
      const bool boost_met =
          boost_err.empty() && (rb.first == rb.second || tol.met(rb.first, rb.second));
      return Rcpp::List::create(
          Rcpp::Named("boost_xs") = boost_xs,
          Rcpp::Named("odelia_xs") = odelia_xs,
          Rcpp::Named("boost_bracket") = Rcpp::NumericVector::create(rb.first, rb.second),
          Rcpp::Named("odelia_bracket") = Rcpp::NumericVector::create(ro.first, ro.second),
          Rcpp::Named("boost_iterations") = static_cast<double>(boost_it),
          Rcpp::Named("odelia_iterations") = static_cast<double>(odelia_it),
          Rcpp::Named("boost_mid") = (rb.first + rb.second) / 2.0,
          Rcpp::Named("odelia_mid") = odelia_mid,
          Rcpp::Named("boost_met") = boost_met,
          Rcpp::Named("boost_err") = boost_err,
          Rcpp::Named("odelia_err") = odelia_err);
    }
  ', env = globalenv())
  .odelia_test_cache$root_find_loaded <- TRUE
}

# Problem, parameter values and bracket, after Alefeld, Potra and Shi (1995),
# with a step (problem 0) beside them.
alefeld_problems <- function() {
  rows <- list(
    list(k = 1, n = 0, a = pi / 2, b = pi),
    list(k = 5, n = 0, a = 0, b = 1.5),
    list(k = 13, n = 0, a = -1, b = 4),
    list(k = 0, n = 0, a = 0, b = 1)
  )
  for (i in 1:10) rows[[length(rows) + 1]] <-
    list(k = 2, n = 0, a = i^2 + 1e-9, b = (i + 1)^2 - 1e-9)
  for (n in c(-40, -100, -200)) rows[[length(rows) + 1]] <-
    list(k = 3, n = n, a = -9, b = 31)
  for (n in c(4, 6, 8, 10, 12)) rows[[length(rows) + 1]] <-
    list(k = 4, n = n, a = 0, b = 5)
  for (n in c(1:5, 20, 40, 60, 80, 100)) rows[[length(rows) + 1]] <-
    list(k = 6, n = n, a = 0, b = 1)
  for (n in c(5, 10, 20)) rows[[length(rows) + 1]] <-
    list(k = 7, n = n, a = 0, b = 1)
  for (n in c(2, 5, 10, 15, 20)) rows[[length(rows) + 1]] <-
    list(k = 8, n = n, a = 0, b = 1)
  for (n in c(1, 2, 4, 5, 8, 15, 20)) rows[[length(rows) + 1]] <-
    list(k = 9, n = n, a = 0, b = 1)
  for (n in c(1, 5, 10, 15, 20)) rows[[length(rows) + 1]] <-
    list(k = 10, n = n, a = 0, b = 1)
  for (n in c(2, 5, 15, 20)) rows[[length(rows) + 1]] <-
    list(k = 11, n = n, a = 0.01, b = 1)
  for (n in c(2, 3, 7, 15, 33)) rows[[length(rows) + 1]] <-
    list(k = 12, n = n, a = 1, b = 100)
  for (n in c(1, 2, 3, 40)) rows[[length(rows) + 1]] <-
    list(k = 14, n = n, a = -1e4, b = pi / 2)
  for (n in c(20, 30, 40, 1000)) rows[[length(rows) + 1]] <-
    list(k = 15, n = n, a = -1e4, b = 1e-4)
  rows
}

# Tolerances as (atol, rtol), and caps from far more than a search needs down
# to few enough that most searches run out.
root_find_settings <- function() {
  expand.grid(tol = c(1e-6, 1e-10, 1e-12, 1e-15), relative_only = c(FALSE, TRUE),
              cap = c(1000, 12, 5))
}

# Every search in the battery under `form`, as the labels of those where odelia
# did not do what Boost did, and how many of the rest converged and ran out.
compare_searches <- function(form) {
  settings <- root_find_settings()
  differ <- character()
  met <- 0L
  ran_out <- 0L
  for (p in alefeld_problems()) {
    for (i in seq_len(nrow(settings))) {
      s <- settings[i, ]
      atol <- if (s$relative_only) 0 else s$tol
      r <- root_find_both(form, p$k, p$n, p$a, p$b, atol, s$tol, s$cap)
      same <- identical(r$odelia_xs, r$boost_xs) && !nzchar(r$boost_err)
      if (form == "bracket") {
        same <- same && identical(r$odelia_bracket, r$boost_bracket) &&
          identical(r$odelia_iterations, r$boost_iterations)
      }
      if (r$boost_met) {
        met <- met + 1L
        same <- same && identical(r$odelia_err, "") &&
          identical(r$odelia_mid, r$boost_mid)
      } else {
        ran_out <- ran_out + 1L
        same <- same && grepl("^RootFindError", r$odelia_err)
      }
      if (!same) {
        differ <- c(differ, sprintf(
          "problem %d (n = %g) on [%g, %g], tol %g%s, cap %d", p$k, p$n, p$a,
          p$b, s$tol, if (s$relative_only) " relative" else "", s$cap))
      }
    }
  }
  list(differ = differ, met = met, ran_out = ran_out)
}

test_that("TOMS748 from the caller's end values is Boost's search", {
  skip_if_not_installed("BH")
  compile_root_find_comparison()
  r <- compare_searches("bracket")
  expect_identical(r$differ, character())
  expect_gt(r$met, 100)
  expect_gt(r$ran_out, 100)
})

test_that("TOMS748 evaluating the ends is Boost's search", {
  skip_if_not_installed("BH")
  compile_root_find_comparison()
  r <- compare_searches("ends")
  expect_identical(r$differ, character())
  expect_gt(r$met, 100)
  expect_gt(r$ran_out, 100)
})

test_that("bisection is Boost's search", {
  skip_if_not_installed("BH")
  compile_root_find_comparison()
  r <- compare_searches("bisect")
  expect_identical(r$differ, character())
  expect_gt(r$met, 100)
  expect_gt(r$ran_out, 100)
})

test_that("ends that bracket no root are refused where Boost refuses them", {
  skip_if_not_installed("BH")
  compile_root_find_comparison()
  for (form in c("bracket", "ends", "bisect")) {
    # x^2 - 0.2 is positive on [1, 2]; the same ends reversed are out of order.
    for (ab in list(c(1, 2), c(2, 1))) {
      r <- root_find_both(form, 4, 2, ab[[1]], ab[[2]], 1e-12, 1e-12, 1000)
      label <- sprintf("%s on [%g, %g]", form, ab[[1]], ab[[2]])
      expect_true(nzchar(r$boost_err), label = label)
      expect_match(r$odelia_err, "^BracketError", label = label)
    }
  }
})
