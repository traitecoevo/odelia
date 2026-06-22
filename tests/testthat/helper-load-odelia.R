.odelia_test_cache <- new.env(parent = emptyenv())
.odelia_test_cache$ode_loaded <- FALSE
.odelia_test_cache$leaf_loaded <- FALSE
.odelia_test_cache$odelia_so <- NA_character_

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

  odelia_so <- NA_character_
  testthat::expect_no_error(
    odelia_so <- odelia_load_dll(local = FALSE, now = TRUE)
  )
  .odelia_test_cache$odelia_so <- odelia_so

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
  #
  # The XAD Tape<T,N> template methods are explicitly instantiated only in
  # src/Tape.cpp, which is compiled into the odelia shared library. The leaf
  # interface is built standalone by sourceCpp, so those symbols must be
  # resolved against the odelia library. On Linux/macOS odelia_load_dll() makes
  # them globally visible (RTLD_GLOBAL) and they resolve at load time, but
  # Windows has no global symbol namespace - DLL imports must be resolved at
  # link time. Linking the sourceCpp build directly against the odelia library
  # via PKG_LIBS (honoured by R CMD SHLIB) works on every platform.
  pkg_cppflags <- paste0("-I", include_dir)
  odelia_so <- .odelia_test_cache$odelia_so
  pkg_libs <- if (is.character(odelia_so) &&
                  length(odelia_so) == 1 &&
                  !is.na(odelia_so) &&
                  nzchar(odelia_so) &&
                  file.exists(odelia_so)) {
    shQuote(normalizePath(odelia_so, winslash = "/", mustWork = FALSE))
  } else {
    Sys.getenv("PKG_LIBS", unset = "")
  }
  withr::local_envvar(
    PKG_CPPFLAGS = pkg_cppflags,
    PKG_LIBS = pkg_libs
  )

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
