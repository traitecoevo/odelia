library(testthat)
invisible(requireNamespace("odelia", quietly = TRUE))
if (!("package:odelia" %in% search())) {
  suppressPackageStartupMessages(
    suppressWarnings(
      library(odelia)
    )
  )
}

.odelia_test_cache <- new.env(parent = emptyenv())
.odelia_test_cache$ode_loaded <- FALSE
.odelia_test_cache$leaf_loaded <- FALSE

is_pkgload_dll <- function() {
  loaded <- getLoadedDLLs()
  if (!("odelia" %in% names(loaded))) {
    return(FALSE)
  }
  dll_path <- tryCatch(loaded[["odelia"]][["path"]], error = function(e) "")
  is.character(dll_path) && length(dll_path) == 1 && grepl("pkgload", dll_path, fixed = TRUE)
}

ensure_ode_interface_loaded <- function(rebuild = FALSE) {
  if (!rebuild && isTRUE(.odelia_test_cache$ode_loaded)) {
    return(invisible(TRUE))
  }

  testthat::expect_no_error(odelia_load_dll(local = FALSE, now = TRUE))

  .odelia_test_cache$ode_loaded <- TRUE
  invisible(TRUE)
}

ensure_leaf_thermal_interfaces <- function(rebuild = FALSE) {
  if (!rebuild && isTRUE(.odelia_test_cache$leaf_loaded)) {
    return(invisible(TRUE))
  }

  ensure_ode_interface_loaded(rebuild = rebuild)

  withr::local_envvar(PKG_CPPFLAGS = paste0("-I", here::here("inst/include")))

  source_cpp_result <- tryCatch(
    {
      Rcpp::sourceCpp(
        here::here("inst/examples/leaf_thermal/src/leaf_thermal_interface.cpp"),
        rebuild = rebuild,
        verbose = FALSE
      )
      NULL
    },
    error = function(e) e
  )

  if (inherits(source_cpp_result, "error")) {
    msg <- conditionMessage(source_cpp_result)
    if (grepl("active_tape_", msg, fixed = TRUE)) {
      testthat::skip("Leaf thermal sourceCpp symbols are unavailable in this load_all session; run installed-package tests for this context.")
    }
    stop(source_cpp_result)
  }

  testthat::expect_silent(
    source(here::here("inst/examples/leaf_thermal/R/leaf_thermal_interface.R"))
  )

  .odelia_test_cache$leaf_loaded <- TRUE
  invisible(TRUE)
}
