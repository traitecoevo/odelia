# The shipped one-file example, inst/examples/consumer_resource/, compiled the
# way the README says to, and checked against the same model as an R function.

ensure_consumer_resource_example <- function() {
  if (isTRUE(.odelia_test_cache$cr_example_loaded)) {
    return(invisible(TRUE))
  }
  ensure_ode_interface_loaded()
  example_cpp <- resolve_test_path(
    "examples/consumer_resource/consumer_resource.cpp",
    "inst/examples/consumer_resource/consumer_resource.cpp")
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
    Rcpp::sourceCpp(example_cpp, verbose = FALSE)
    NULL
  }, error = function(e) e)
  if (inherits(res, "error")) {
    if (grepl("active_tape_", conditionMessage(res), fixed = TRUE)) {
      testthat::skip("consumer-resource example symbols unavailable in this load_all session.")
    }
    stop(res)
  }
  .odelia_test_cache$cr_example_loaded <- TRUE
  invisible(TRUE)
}

cr_example_rhs <- function(t, y, p) {
  uptake <- p[["a"]] * y[2] / (p[["h"]] + y[2])
  list(c(y[1] * (uptake - p[["m"]]),
         p[["r"]] * y[2] * (1 - y[2] / p[["K"]]) - p[["q"]] * y[1] * uptake))
}
cr_example_pars <- c(a = 1, h = 0.5, q = 0.5, m = 0.4, r = 1, K = 1)
cr_example_y0 <- c(N = 0.5, R = 0.5)

testthat::test_that("the one-file example solves the model its R twin solves", {
  ensure_consumer_resource_example()
  times <- seq(0, 40, by = 2)
  for (method in c("dopri", "rodas")) {
    out <- cr_solve(cr_solver(cr_example_pars, cr_example_y0, 0, method, 1e-10), times)
    ref <- ode_solve(cr_example_rhs, cr_example_y0, times, parms = cr_example_pars,
                     method = method, rtol = 1e-10, atol = 1e-10, autonomous = TRUE)
    expect_identical(colnames(out), c("time", "N", "R"))
    expect_equal(out[, "time"], times)
    expect_lt(max(abs(out[, 2:3] - ref[, 2:3])), 1e-8)
  }
})

testthat::test_that("the one-file example's sensitivities are central differences of its solves", {
  ensure_consumer_resource_example()
  times <- seq(0, 40, by = 10)
  for (method in c("dopri", "rodas")) {
    sens <- cr_sensitivities(cr_solver(cr_example_pars, cr_example_y0, 0, method, 1e-8),
                             times, 0:5, 0:1)
    expect_identical(dim(sens$dy), c(length(times), 2L, 8L))
    inputs <- c(cr_example_pars, cr_example_y0)
    for (j in seq_along(inputs)) {
      h <- 1e-5 * max(1, abs(inputs[[j]]))
      up <- inputs; up[[j]] <- up[[j]] + h
      down <- inputs; down[[j]] <- down[[j]] - h
      fd <- (cr_solve(cr_solver(up[1:6], up[7:8], 0, method, 1e-8), times)[, 2:3] -
             cr_solve(cr_solver(down[1:6], down[7:8], 0, method, 1e-8), times)[, 2:3]) / (2 * h)
      expect_lt(max(abs(sens$dy[, , j] - fd) / pmax(1, abs(fd))), 1e-5)
    }
  }
})
