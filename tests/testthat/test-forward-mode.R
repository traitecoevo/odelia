# Forward-mode sensitivities: a run's steps walked again carrying one direction
# per input (sweep_tangent), from R. Refereed against the reverse sweep where
# both apply, and against central differences of the same pinned steps under
# RODAS, which cannot be swept.

ensure_cr_runner <- function() {
  if (isTRUE(.odelia_test_cache$cr_loaded)) {
    return(invisible(TRUE))
  }
  ensure_ode_interface_loaded()
  runner_cpp <- resolve_test_path(
    "tests/testthat/consumer_resource_runner.cpp",
    "tests/testthat/consumer_resource_runner.cpp")
  odelia_so <- .odelia_test_cache$odelia_so
  pkg_libs <- if (is.character(odelia_so) && length(odelia_so) == 1 &&
                  !is.na(odelia_so) && nzchar(odelia_so) &&
                  file.exists(odelia_so)) {
    shQuote(normalizePath(odelia_so, winslash = "/", mustWork = FALSE))
  } else {
    Sys.getenv("PKG_LIBS", unset = "")
  }
  withr::local_envvar(PKG_CPPFLAGS = odelia_cppflags(), PKG_LIBS = pkg_libs)
  res <- tryCatch({
    Rcpp::sourceCpp(runner_cpp, verbose = FALSE)
    NULL
  }, error = function(e) e)
  if (inherits(res, "error")) {
    if (grepl("active_tape_", conditionMessage(res), fixed = TRUE)) {
      testthat::skip("consumer-resource runner symbols unavailable in this load_all session.")
    }
    stop(res)
  }
  .odelia_test_cache$cr_loaded <- TRUE
  invisible(TRUE)
}

# The largest gap between two arrays, relative to each entry's size or 1.
largest_gap <- function(a, b) {
  max(abs(a - b) / pmax(1, abs(b)))
}

# a, h, q, m, r, K. A half-saturation of 1e-3 makes the resource relax a
# thousand times faster than the consumer: stiff.
cr_pars <- c(1, 1e-3, 0.5, 0.2, 1, 2)
cr_y0 <- c(0.5, 1)

testthat::test_that("forward and reverse give one gradient on a stiff system under explicit steppers", {
  ensure_cr_runner()
  times <- seq(0, 5, by = 0.5)
  for (method in c("rkck", "dopri")) {
    schedule <- cr_schedule(cr_pars, cr_y0, times, method, 1e-8)
    obs <- match(times[-1], schedule)
    target <- cr_pinned(cr_pars, cr_y0, schedule, method)[obs, ] * 1.1
    s <- cr_solver(cr_pars, cr_y0, method, 1e-8)
    # At the parameters the schedule was chosen at: an explicit step sized for
    # them sits near its stability limit on this system, and pinned steps at
    # other parameters can cross it.
    fwd <- cr_loss_and_gradient(s, schedule, target, obs, cr_y0, cr_pars, "forward")
    rev <- cr_loss_and_gradient(s, schedule, target, obs, cr_y0, cr_pars, "reverse")
    expect_equal(fwd$loss, rev$loss, tolerance = 1e-12)
    expect_length(fwd$gradient, 8)
    expect_lt(largest_gap(fwd$gradient, rev$gradient), 1e-10)
  }
})

testthat::test_that("RODAS sensitivities are the central difference of the same pinned RODAS steps", {
  ensure_cr_runner()
  times <- seq(0, 5, by = 0.5)
  s <- cr_solver(cr_pars, cr_y0, "rodas", 1e-8)
  sens <- cr_sensitivities(s, times, 0:5, 0:1)
  expect_identical(dim(sens$dy), c(length(times), 2L, 8L))

  schedule <- cr_schedule(cr_pars, cr_y0, times, "rodas", 1e-8)
  rows <- match(times, schedule)
  expect_false(anyNA(rows))
  inputs <- c(cr_pars, cr_y0)
  walk <- function(x) cr_pinned(x[1:6], x[7:8], schedule, "rodas")[rows, ]
  for (j in seq_along(inputs)) {
    h <- 1e-6 * max(1, abs(inputs[[j]]))
    up <- inputs; up[[j]] <- up[[j]] + h
    down <- inputs; down[[j]] <- down[[j]] - h
    fd <- (walk(up) - walk(down)) / (2 * h)
    expect_lt(largest_gap(sens$dy[, , j], fd), 1e-6)
  }
  expect_lt(largest_gap(sens$y, walk(inputs)), 1e-9)
})

testthat::test_that("RODAS sensitivities converge on the explicit stepper's as tolerances tighten", {
  ensure_cr_runner()
  times <- c(0, 2)
  reference <- cr_sensitivities(cr_solver(cr_pars, cr_y0, "dopri", 1e-11),
                                times, 0:5, 0:1)$dy[2, , ]
  gaps <- vapply(c(1e-4, 1e-6, 1e-8), function(tol) {
    largest_gap(cr_sensitivities(cr_solver(cr_pars, cr_y0, "rodas", tol),
                                 times, 0:5, 0:1)$dy[2, , ], reference)
  }, numeric(1))
  expect_true(all(diff(gaps) < 0))
  expect_lt(gaps[[3]], 1e-5)
})

testthat::test_that("how many lanes a pass carries changes no derivative", {
  ensure_cr_runner()
  times <- seq(0, 2, by = 0.5)
  for (method in c("dopri", "rodas")) {
    wide <- cr_sensitivities(cr_solver(cr_pars, cr_y0, method, 1e-8), times, 0:5, 0:1)
    narrow <- cr_sensitivities_narrow(cr_solver(cr_pars, cr_y0, method, 1e-8),
                                      times, 0:5, 0:1)
    expect_lt(largest_gap(narrow$dy, wide$dy), 1e-13)
  }
})

lorenz_solver <- function(method = "rkck") {
  lz <- LorenzSystem$new(10, 28, 8 / 3)
  lz$set_initial_state(c(1, 1, 1), 0)
  Lorenz_Solver$new(lz$ptr, odelia:::OdeControl$new()$ptr, method)
}

testthat::test_that("$sensitivities() names its answer, selects inputs, and leaves the solver where it was", {
  s <- lorenz_solver("dopri")
  times <- c(0, 0.25, 0.5)
  out <- s$sensitivities(times)
  expect_identical(out$time, times)
  expect_identical(dim(out$y), c(3L, 3L))
  expect_identical(colnames(out$y), c("x", "y", "z"))
  expect_identical(dim(out$dy), c(3L, 3L, 6L))
  expect_identical(dimnames(out$dy)[[3]], c("sigma", "R", "b", "x0", "y0", "z0"))
  expect_equal(unname(out$dy[1, , 4:6]), diag(3))
  expect_identical(s$time(), 0)

  plain <- lorenz_solver("dopri")
  plain$advance_adaptive(times)
  expect_lt(largest_gap(out$y[3, ], plain$state()), 1e-9)

  some <- s$sensitivities(times, parameters = "R", initial = 3)
  expect_identical(dimnames(some$dy)[[3]], c("R", "z0"))
  expect_identical(some$dy[, , "R"], out$dy[, , "R"])
  expect_identical(some$dy[, , "z0"], out$dy[, , "z0"])
  none <- s$sensitivities(times, parameters = FALSE, initial = FALSE)
  expect_identical(dim(none$dy), c(3L, 3L, 0L))

  expect_error(s$sensitivities(times, parameters = "rho"), "names not among")
  expect_error(s$sensitivities(times, initial = 4), "must lie in 1..3")
  expect_error(s$sensitivities(c(1, 2)), "current time")
})

testthat::test_that("$sensitivities() equals one reverse sweep per output", {
  times <- c(0, 0.5, 1)
  fwd <- lorenz_solver("rkck")$sensitivities(times)
  # The reverse gradient of a loss that observes one component at one time is
  # twice the difference times that component's row of dy.
  reference <- lorenz_solver("rkck")
  reference$advance_adaptive(times)
  schedule <- reference$times()
  obs <- match(times[3], schedule)
  fitter <- lorenz_solver("rkck")
  fitter$set_target(schedule, matrix(fwd$y[3, ] - c(0.5, 0, 0), 1), obs)
  rev <- fitter$loss_and_gradient(initial_state = c(1, 1, 1),
                                  params = c(10, 28, 8 / 3), mode = "reverse")
  expect_lt(largest_gap(rev$gradient, 2 * 0.5 * fwd$dy[3, 1, ]), 1e-10)
})

testthat::test_that("$loss_and_gradient() is forward under rodas and refuses a reverse sweep there", {
  ref <- lorenz_solver("rkck")
  ref$advance_adaptive(seq(0, 1, by = 0.1))
  schedule <- ref$times()
  hist <- ref$history()
  target <- as.matrix(hist[, c("x", "y", "z")])
  obs <- match(hist$time, schedule)
  p <- c(11, 27, 2.5)

  fitter <- function(method) {
    s <- lorenz_solver(method)
    s$set_target(schedule, target, obs)
    s
  }
  for (method in c("rkck", "dopri")) {
    s <- fitter(method)
    fwd <- s$loss_and_gradient(params = p, mode = "forward")
    expect_identical(s$loss_and_gradient(params = p), s$loss_and_gradient(params = p, mode = "reverse"))
    expect_lt(largest_gap(fwd$gradient, s$loss_and_gradient(params = p)$gradient), 1e-10)
  }

  s <- fitter("rodas")
  res <- s$loss_and_gradient(params = p)
  expect_identical(res, s$loss_and_gradient(params = p, mode = "forward"))
  fd <- vapply(seq_along(p), function(i) {
    h <- 1e-6 * max(1, abs(p[[i]]))
    up <- p; dn <- p
    up[[i]] <- up[[i]] + h
    dn[[i]] <- dn[[i]] - h
    (s$loss_and_gradient(params = up)$loss - s$loss_and_gradient(params = dn)$loss) / (2 * h)
  }, numeric(1))
  expect_lt(largest_gap(res$gradient, fd), 1e-6)
  expect_error(s$loss_and_gradient(params = p, mode = "reverse"),
               "cannot be swept in reverse")
})
