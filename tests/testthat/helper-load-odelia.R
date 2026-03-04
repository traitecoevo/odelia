library(testthat)
suppressPackageStartupMessages(
  suppressWarnings(
    library(odelia)
  )
)

.odelia_test_cache <- new.env(parent = emptyenv())
.odelia_test_cache$ode_loaded <- FALSE
.odelia_test_cache$leaf_loaded <- FALSE

ensure_ode_interface_loaded <- function(rebuild = FALSE) {
  if (!rebuild && isTRUE(.odelia_test_cache$ode_loaded)) {
    return(invisible(TRUE))
  }

  withr::local_envvar(PKG_CPPFLAGS = paste0("-I", here::here("inst/include")))

  odelia_so <- file.path(
    system.file("libs", .Platform$r_arch, package = "odelia"),
    paste0("odelia", .Platform$dynlib.ext)
  )
  testthat::expect_true(file.exists(odelia_so))
  testthat::expect_silent(
    suppressWarnings(try(dyn.load(odelia_so, local = FALSE, now = TRUE), silent = TRUE))
  )

  testthat::expect_no_error(
    Rcpp::sourceCpp(here::here("src/ode_interface.cpp"), rebuild = rebuild, verbose = FALSE)
  )

  .odelia_test_cache$ode_loaded <- TRUE
  invisible(TRUE)
}

ensure_leaf_thermal_interfaces <- function(rebuild = FALSE) {
  if (!rebuild && isTRUE(.odelia_test_cache$leaf_loaded)) {
    return(invisible(TRUE))
  }

  ensure_ode_interface_loaded(rebuild = rebuild)

  testthat::expect_no_error(
    Rcpp::sourceCpp(
      here::here("inst/examples/leaf_thermal/src/leaf_thermal_interface.cpp"),
      rebuild = rebuild,
      verbose = FALSE
    )
  )
  testthat::expect_silent(
    source(here::here("inst/examples/leaf_thermal/R/leaf_thermal_interface.R"))
  )

  .odelia_test_cache$leaf_loaded <- TRUE
  invisible(TRUE)
}
