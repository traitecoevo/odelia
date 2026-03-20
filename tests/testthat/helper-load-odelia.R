.odelia_test_cache <- new.env(parent = emptyenv())
.odelia_test_cache$ode_loaded <- FALSE
.odelia_test_cache$leaf_loaded <- FALSE

resolve_test_path <- function(installed_rel, source_rel) {
  pkg_path <- system.file(package = "odelia")
  if (nzchar(pkg_path)) {
    installed_path <- file.path(pkg_path, installed_rel)
    if (file.exists(installed_path)) {
      return(installed_path)
    }
  }

  source_path <- here::here(source_rel)
  if (file.exists(source_path)) {
    return(source_path)
  }

  stop(
    sprintf(
      "Could not locate required test file: installed='%s', source='%s'",
      installed_rel,
      source_rel
    )
  )
}

# Detect whether the currently loaded DLL came from pkgload/load_all,
# which can behave differently from an installed package binary.
is_pkgload_dll <- function() {
  loaded <- getLoadedDLLs()
  # If the package DLL isn't loaded, we can't be in a pkgload session.
  if (!("odelia" %in% names(loaded))) {
    return(FALSE)
  }
  # Check if the loaded DLL path contains "pkgload", which is a common indicator of a pkgload session.
  dll_path <- tryCatch(loaded[["odelia"]][["path"]], error = function(e) "")
  # Ensure the path is a single string and contains "pkgload" to confirm we're in a pkgload session.
  is.character(dll_path) && length(dll_path) == 1 && grepl("pkgload", dll_path, fixed = TRUE)
}

# Ensure the package-level ODE interface symbols are loaded exactly once
# per test session unless an explicit rebuild is requested.
ensure_ode_interface_loaded <- function(rebuild = FALSE) {
  if (!rebuild && isTRUE(.odelia_test_cache$ode_loaded)) {
    return(invisible(TRUE))
  }

  testthat::expect_no_error(odelia_load_dll(local = FALSE, now = TRUE))

  .odelia_test_cache$ode_loaded <- TRUE
  invisible(TRUE)
}

# Compile and source the leaf thermal test interface on demand, while
# gracefully skipping known load_all symbol limitations for AD symbols.
ensure_leaf_thermal_interfaces <- function(rebuild = FALSE) {
  if (!rebuild && isTRUE(.odelia_test_cache$leaf_loaded)) {
    return(invisible(TRUE))
  }

  ensure_ode_interface_loaded(rebuild = rebuild)

  include_dir <- dirname(dirname(resolve_test_path("include/odelia/ode_solver.hpp", "inst/include/odelia/ode_solver.hpp")))
  leaf_cpp <- resolve_test_path(
    "examples/leaf_thermal/src/leaf_thermal_interface.cpp",
    "inst/examples/leaf_thermal/src/leaf_thermal_interface.cpp"
  )
  leaf_r <- resolve_test_path(
    "examples/leaf_thermal/R/leaf_thermal_interface.R",
    "inst/examples/leaf_thermal/R/leaf_thermal_interface.R"
  )

  # Provide include paths needed by sourceCpp for package headers.
  withr::local_envvar(PKG_CPPFLAGS = paste0("-I", include_dir))

  source_cpp_result <- tryCatch(
    {
      Rcpp::sourceCpp(
        leaf_cpp,
        rebuild = rebuild,
        verbose = FALSE
      )
      NULL
    },
    error = function(e) e
  )

  if (inherits(source_cpp_result, "error")) {
    msg <- conditionMessage(source_cpp_result)
    # In dev load_all contexts, active_tape_ symbols may be unavailable.
    # Skip here to avoid false negatives; installed-package tests cover it.
    if (grepl("active_tape_", msg, fixed = TRUE)) {
      testthat::skip("Leaf thermal sourceCpp symbols are unavailable in this load_all session; run installed-package tests for this context.")
    }
    stop(source_cpp_result)
  }

  testthat::expect_silent(
    source(leaf_r)
  )

  .odelia_test_cache$leaf_loaded <- TRUE
  invisible(TRUE)
}
