# Tests for odelia::ode::vector_jacobian_product: one block recorded and swept
# once, writing input adjoints into a caller-owned buffer.
#
# The XAD adjoint Tape<T,N> methods are instantiated only in src/Tape.cpp (in the
# odelia shared library), so this sourceCpp build links against it via PKG_LIBS,
# exactly like the leaf-thermal AD interface (see helper-load-odelia.R).

compile_vjp_interface <- function() {
  ensure_ode_interface_loaded(rebuild = FALSE)
  include_dir <- odelia_include_dir()
  odelia_so <- .odelia_test_cache$odelia_so
  pkg_libs <- if (is.character(odelia_so) && length(odelia_so) == 1 &&
                  !is.na(odelia_so) && nzchar(odelia_so) && file.exists(odelia_so)) {
    shQuote(normalizePath(odelia_so, winslash = "/", mustWork = FALSE))
  } else {
    Sys.getenv("PKG_LIBS", unset = "")
  }
  withr::local_envvar(PKG_CPPFLAGS = odelia_cppflags(include_dir),
                      PKG_LIBS = pkg_libs)

  res <- tryCatch({
    Rcpp::sourceCpp(code = '
      // [[Rcpp::plugins(cpp20)]]
      #include <Rcpp.h>
      #include <cstddef>
      #include <vector>
      #include <XAD/XAD.hpp>
      #include <odelia/adjoint.hpp>

      // A four-input, three-output block, small enough to difference centrally.
      template <class S>
      static void small_block(const std::vector<S>& x, std::vector<S>& y) {
        y[0] = sin(x[0]) * x[1] + x[2] * x[2] * x[3];
        y[1] = exp(0.5 * x[1]) - x[0] * x[3];
        y[2] = x[0] * x[1] * x[2] + log(1.0 + x[3] * x[3]);
      }

      // Fixed arithmetic on the first three entries; the rest are registered inputs
      // the block never reads, so the recording does not grow with their number.
      template <class S>
      static void padded_block(const std::vector<S>& x, std::vector<S>& y) {
        for (std::size_t i = 0; i < y.size(); ++i) {
          S v = sin(x[0]) * x[1] + x[2] * x[2] + S(double(i)) * x[1];
          y[i] = v;
        }
      }

      // [[Rcpp::export]]
      Rcpp::NumericVector small_block_value(std::vector<double> x) {
        std::vector<double> y(3);
        small_block(x, y);
        return Rcpp::wrap(y);
      }

      // [[Rcpp::export]]
      Rcpp::NumericVector small_block_vjp(std::vector<double> x,
                                          std::vector<double> output_adjoints) {
        xad::adj<double>::tape_type tape(false);
        const odelia::ode::adjoint_rows seeds =
            odelia::ode::adjoint_rows::one_row(output_adjoints);
        odelia::ode::adjoint_rows input_adjoints;
        odelia::ode::vector_jacobian_product(
            tape, x, seeds,
            [](const auto& xa, auto& ya) { small_block(xa, ya); }, input_adjoints);
        return Rcpp::wrap(input_adjoints.to_rows()[0]);
      }

      // [[Rcpp::export]]
      double padded_block_recording_size(int n_inputs,
                                         std::vector<double> output_adjoints) {
        std::vector<double> x(n_inputs, 0.5);
        xad::adj<double>::tape_type tape(false);
        const odelia::ode::adjoint_rows seeds =
            odelia::ode::adjoint_rows::one_row(output_adjoints);
        odelia::ode::adjoint_rows input_adjoints;
        return double(odelia::ode::vector_jacobian_product(
            tape, x, seeds,
            [](const auto& xa, auto& ya) { padded_block(xa, ya); }, input_adjoints));
      }

      // Repeated products on one caller-owned tape. Returns the adjoints of every
      // call laid end to end, then the tape memory after each call: a tape that is
      // not reset between calls either carries adjoints over from the call before it
      // or grows its recording, and both show up here.
      // [[Rcpp::export]]
      Rcpp::List small_block_vjp_reused_tape(std::vector<double> x,
                                             std::vector<double> output_adjoints,
                                             int n_calls) {
        xad::adj<double>::tape_type tape(false);
        const odelia::ode::adjoint_rows seeds =
            odelia::ode::adjoint_rows::one_row(output_adjoints);
        odelia::ode::adjoint_rows input_adjoints;
        std::vector<double> all_adjoints, sizes;
        for (int k = 0; k < n_calls; ++k) {
          sizes.push_back(double(odelia::ode::vector_jacobian_product(
              tape, x, seeds,
              [](const auto& xa, auto& ya) { small_block(xa, ya); }, input_adjoints)));
          for (double v : input_adjoints[0]) all_adjoints.push_back(v);
        }
        return Rcpp::List::create(Rcpp::Named("adjoints") = Rcpp::wrap(all_adjoints),
                                  Rcpp::Named("sizes") = Rcpp::wrap(sizes));
      }

      // [[Rcpp::export]]
      Rcpp::NumericVector small_block_vjp_under_active_tape(
          std::vector<double> x, std::vector<double> output_adjoints) {
        xad::adj<double>::tape_type tape;
        return small_block_vjp(x, output_adjoints);
      }

      // A caller-owned tape handed in while a second tape is the active one. The
      // product must record onto neither.
      // [[Rcpp::export]]
      Rcpp::NumericVector small_block_vjp_owned_tape_under_another(
          std::vector<double> x, std::vector<double> output_adjoints) {
        xad::adj<double>::tape_type owned(false);
        xad::adj<double>::tape_type other;  // constructed active
        const odelia::ode::adjoint_rows seeds =
            odelia::ode::adjoint_rows::one_row(output_adjoints);
        odelia::ode::adjoint_rows input_adjoints;
        odelia::ode::vector_jacobian_product(
            owned, x, seeds,
            [](const auto& xa, auto& ya) { small_block(xa, ya); }, input_adjoints);
        return Rcpp::wrap(input_adjoints.to_rows()[0]);
      }', verbose = FALSE)
    NULL
  }, error = function(e) e)

  if (inherits(res, "error")) {
    msg <- conditionMessage(res)
    if (grepl("active_tape_", msg, fixed = TRUE)) {
      testthat::skip("vector_jacobian_product sourceCpp symbols are unavailable in this load_all session; run installed-package tests for this context.")
    }
    stop(res)
  }
  invisible(TRUE)
}

testthat::test_that("vector_jacobian_product reproduces a central finite difference", {
  compile_vjp_interface()

  x <- c(0.7, -1.3, 0.45, 1.1)
  w <- c(1.0, -2.0, 0.5)

  # The product is t(J) %*% w, so its j-th entry is the derivative of w . f(x)
  # along input j -- which a central difference of the same block gives directly.
  h  <- 1e-5
  fd <- vapply(seq_along(x), function(j) {
    xp <- x; xp[j] <- xp[j] + h
    xm <- x; xm[j] <- xm[j] - h
    sum(w * (small_block_value(xp) - small_block_value(xm))) / (2 * h)
  }, numeric(1))

  expect_equal(small_block_vjp(x, w), fd, tolerance = 1e-8)
})

testthat::test_that("the recording size is a property of the block", {
  compile_vjp_interface()

  # One recording serves every output adjoint, so sweeping the same three-output
  # block with one seeded adjoint and with three gives the same recording size.
  # A product implemented as a loop of gradients would record once per output.
  expect_equal(padded_block_recording_size(10, c(1, 0, 0)),
               padded_block_recording_size(10, c(1, 2, 3)))

  # Inputs the block does not read add registered slots, not recorded operations,
  # so a tenfold input count must not multiply the recording size tenfold.
  small <- padded_block_recording_size(10, c(1, 2, 3))
  large <- padded_block_recording_size(100, c(1, 2, 3))
  expect_lt(large, 10 * small)
})

testthat::test_that("a reused tape gives each call the adjoints a single call gives", {
  compile_vjp_interface()

  x <- c(0.7, -1.3, 0.45, 1.1)
  w <- c(1.0, -2.0, 0.5)
  once <- small_block_vjp(x, w)

  n <- 4L
  reused <- small_block_vjp_reused_tape(x, w, n)

  # Every call on the shared tape sees an empty recording, so it reports the same
  # adjoints as a call that owned its own tape.
  expect_equal(reused$adjoints, rep(once, n))

  # And the recording is the same size on every call. Without the reset the tape
  # keeps the previous calls' derivative slots and this grows without bound.
  expect_equal(reused$sizes, rep(reused$sizes[1], n))
})

testthat::test_that("vector_jacobian_product stops when a tape is already active", {
  compile_vjp_interface()

  expect_error(small_block_vjp_under_active_tape(c(0.7, -1.3, 0.45, 1.1),
                                                 c(1.0, -2.0, 0.5)),
               "tape is already active")

  # Handing a tape in does not license recording onto whichever tape happens to be
  # active: the one handed in must be the active one.
  expect_error(small_block_vjp_owned_tape_under_another(c(0.7, -1.3, 0.45, 1.1),
                                                        c(1.0, -2.0, 0.5)),
               "tape is already active")
})
