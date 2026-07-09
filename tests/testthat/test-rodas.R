# Tests for the implicit RODAS4(3) Rosenbrock stepper (method = "rodas").

# Compile the standalone Van der Pol runner on demand, wiring include/link flags
# the same way the leaf-thermal helper does (see helper-load-odelia.R).
ensure_vdp_runner <- function() {
  if (isTRUE(.odelia_test_cache$vdp_loaded)) {
    return(invisible(TRUE))
  }
  ensure_ode_interface_loaded()

  include_dir <- dirname(dirname(resolve_test_path(
    "include/odelia/ode_solver.hpp", "inst/include/odelia/ode_solver.hpp")))
  runner_cpp <- resolve_test_path(
    "tests/testthat/vanderpol_runner.cpp",
    "tests/testthat/vanderpol_runner.cpp")

  odelia_so <- .odelia_test_cache$odelia_so
  pkg_libs <- if (is.character(odelia_so) && length(odelia_so) == 1 &&
                  !is.na(odelia_so) && nzchar(odelia_so) &&
                  file.exists(odelia_so)) {
    shQuote(normalizePath(odelia_so, winslash = "/", mustWork = FALSE))
  } else {
    Sys.getenv("PKG_LIBS", unset = "")
  }
  withr::local_envvar(
    PKG_CPPFLAGS = paste0("-I", include_dir),
    PKG_LIBS = pkg_libs
  )

  res <- tryCatch({
    Rcpp::sourceCpp(runner_cpp, verbose = FALSE)
    NULL
  }, error = function(e) e)
  if (inherits(res, "error")) {
    if (grepl("active_tape_", conditionMessage(res), fixed = TRUE)) {
      testthat::skip("VdP runner symbols unavailable in this load_all session.")
    }
    stop(res)
  }

  .odelia_test_cache$vdp_loaded <- TRUE
  invisible(TRUE)
}

lorenz_r_rhs <- function(t, y, p) {
  list(c(p[["sigma"]] * (y[[2]] - y[[1]]),
         p[["R"]] * y[[1]] - y[[2]] - y[[1]] * y[[3]],
         -p[["b"]] * y[[3]] + y[[1]] * y[[2]]))
}

testthat::test_that("RODAS and RKCK agree on the (non-stiff) Lorenz system", {
  ensure_ode_interface_loaded()
  pars <- c(sigma = 10.0, R = 28.0, b = 8.0 / 3.0)
  times <- seq(0, 2, by = 0.05)

  ctrl <- odelia:::OdeControl$new()
  ctrl$set_tol_rel(1e-10)
  ctrl$set_tol_abs(1e-10)

  run <- function(method) {
    lz <- odelia:::LorenzSystem$new(pars[["sigma"]], pars[["R"]], pars[["b"]])
    lz$set_state(c(1, 1, 1), 0.0)
    runner <- odelia:::Lorenz_Solver$new(lz$ptr, ctrl$ptr, active = FALSE,
                                         method = method)
    runner$advance_adaptive(times)
    runner$history()
  }

  out_rkck <- run("rkck")
  out_rodas <- run("rodas")

  expect_equal(nrow(out_rodas), length(times))
  expect_true(all(is.finite(out_rodas$x)))
  # Both are accurate integrators of the same smooth trajectory over a short,
  # non-chaotic-divergence horizon.
  expect_equal(out_rodas$x, out_rkck$x, tolerance = 1e-5)
  expect_equal(out_rodas$y, out_rkck$y, tolerance = 1e-5)
  expect_equal(out_rodas$z, out_rkck$z, tolerance = 1e-5)
})

testthat::test_that("RODAS on Lorenz matches deSolve", {
  skip_if_not_installed("deSolve")
  ensure_ode_interface_loaded()
  pars <- c(sigma = 10.0, R = 28.0, b = 8.0 / 3.0)
  times <- seq(0, 2, by = 0.05)

  ctrl <- odelia:::OdeControl$new()
  ctrl$set_tol_rel(1e-10)
  ctrl$set_tol_abs(1e-10)
  lz <- odelia:::LorenzSystem$new(pars[["sigma"]], pars[["R"]], pars[["b"]])
  lz$set_state(c(1, 1, 1), 0.0)
  runner <- odelia:::Lorenz_Solver$new(lz$ptr, ctrl$ptr, active = FALSE,
                                       method = "rodas")
  runner$advance_adaptive(times)
  out <- runner$history()

  ref <- deSolve::ode(y = c(1, 1, 1), times = times, func = lorenz_r_rhs,
                      parms = pars, method = "radau",
                      rtol = 1e-10, atol = 1e-10)

  expect_equal(out$x, as.numeric(ref[, 2]), tolerance = 1e-5)
  expect_equal(out$y, as.numeric(ref[, 3]), tolerance = 1e-5)
  expect_equal(out$z, as.numeric(ref[, 4]), tolerance = 1e-5)
})

testthat::test_that("RODAS solves stiff Van der Pol and matches deSolve", {
  skip_if_not_installed("deSolve")
  ensure_vdp_runner()

  eps <- 1e-2
  y0 <- c(2, 0)
  times <- seq(0, 2, by = 0.1)

  res <- vdp_run(eps, times, y0, tol_rel = 1e-8, tol_abs = 1e-8, method = "rodas")
  expect_equal(nrow(res$state), length(times))
  expect_true(all(is.finite(res$state)))

  vdp_r <- function(t, y, eps) {
    list(c(y[[2]], ((1 - y[[1]]^2) * y[[2]] - y[[1]]) / eps))
  }
  ref <- deSolve::ode(y = y0, times = times, func = vdp_r, parms = eps,
                      method = "radau", rtol = 1e-8, atol = 1e-8)

  expect_equal(res$state[, 1], as.numeric(ref[, 2]), tolerance = 1e-4)
  expect_equal(res$state[, 2], as.numeric(ref[, 3]), tolerance = 1e-3)
})

testthat::test_that("RODAS takes far fewer steps than RKCK on a stiff problem", {
  ensure_vdp_runner()
  eps <- 1e-4
  y0 <- c(2, 0)
  times <- seq(0, 2, by = 0.2)

  rodas <- vdp_run(eps, times, y0, 1e-6, 1e-6, "rodas")
  rkck  <- vdp_run(eps, times, y0, 1e-6, 1e-6, "rkck")

  expect_true(all(is.finite(rodas$state)))
  # The whole point of an implicit method: stability without step collapse. At
  # eps = 1e-4 the explicit RKCK step is stability-limited (~1e4 steps) while
  # RODAS stays accuracy-limited (~1e3); the gap widens further as eps shrinks.
  expect_lt(rodas$n_steps, rkck$n_steps / 5)
})

testthat::test_that("RODAS is rejected for AD/active solvers (for now)", {
  ensure_ode_interface_loaded()
  ctrl <- odelia:::OdeControl$new()
  lz <- odelia:::LorenzSystem$new(10, 28, 8 / 3)
  lz$set_initial_state(c(1, 1, 1), 0)
  runner <- odelia:::Lorenz_Solver$new(lz$ptr, ctrl$ptr, active = TRUE,
                                       method = "rodas")
  expect_error(runner$advance_fixed(seq(0, 1, by = 0.1)),
               regexp = "rodas")
})
