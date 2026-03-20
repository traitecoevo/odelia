#' Load odelia shared library with configurable symbol visibility
#'
#' Locates the installed 
#' `odelia` shared library and calls [base::dyn.load()] on it.
#' This is useful in `Rcpp::sourceCpp()` workflows where symbols from the
#' package shared object must be available to a temporary shared object.
#'
#' @param local Passed to [base::dyn.load()]. Default `FALSE` to expose symbols
#'   globally for downstream dynamic linking.
#' @param now Passed to [base::dyn.load()]. Default `TRUE`.
#'
#' @return Invisibly returns the resolved path to the `odelia` shared library.
#' @export
odelia_load_dll <- function(local = FALSE, now = TRUE) {
  
  loaded <- getLoadedDLLs()
  is_loaded <- "odelia" %in% names(loaded)
  loaded_path <- NA_character_
  if (is_loaded) {
    loaded_path <- tryCatch(loaded[["odelia"]][["path"]], error = function(e) NA_character_)
    if (is.null(loaded_path) || !is.character(loaded_path) || length(loaded_path) != 1) {
      loaded_path <- NA_character_
    }
  }

  dll_name <- paste0("odelia", .Platform$dynlib.ext)
  lib_dirs <- unique(Filter(
    nzchar,
    c(
      system.file("libs", .Platform$r_arch, package = "odelia"),
      system.file("libs", package = "odelia")
    )
  ))

  candidates <- file.path(lib_dirs, dll_name)
  odelia_so <- candidates[file.exists(candidates)][1]

  if ((is.na(odelia_so) || !nzchar(odelia_so)) && !is.na(loaded_path) && nzchar(loaded_path) && file.exists(loaded_path)) {
    odelia_so <- loaded_path
  }

  if (is.na(odelia_so) || !nzchar(odelia_so) || !file.exists(odelia_so)) {
    stop("Could not find odelia shared library; run make install (or reinstall odelia) and try again.")
  }

  if (!is_loaded) {
    dyn.load(odelia_so, local = local, now = now)
  } else if (!isTRUE(local) && !isTRUE(getOption("odelia.dll.global", FALSE))) {
    dyn.load(odelia_so, local = FALSE, now = now)
    options(odelia.dll.global = TRUE)
  }

  invisible(odelia_so)
}
