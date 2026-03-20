#' Namespace exports for low-level wrappers
#'
#' @description
#' Ensure selected low-level Rcpp wrapper functions remain exported.
#' @name low_level_exports
#' @noRd
#'
#' @keywords internal
#' @rawNamespace export(Drivers_new)
#' @rawNamespace export(Drivers_set_constant)
#' @rawNamespace export(Drivers_set_variable)
#' @rawNamespace export(Drivers_set_extrapolate)
#' @rawNamespace export(Drivers_evaluate)
#' @rawNamespace export(Drivers_evaluate_range)
#' @rawNamespace export(Drivers_get_names)
#' @rawNamespace export(Drivers_clear)
#' @rawNamespace export(OdeControl_new)
#' @rawNamespace export(OdeControl_get_controls)
#' @rawNamespace export(OdeControl_set_controls)
#' @rawNamespace export(OdeControl_set_tol_abs)
#' @rawNamespace export(OdeControl_set_tol_rel)
#' @rawNamespace export(OdeControl_set_a_y)
#' @rawNamespace export(OdeControl_set_a_dydt)
#' @rawNamespace export(OdeControl_set_step_size_min)
#' @rawNamespace export(OdeControl_set_step_size_max)
#' @rawNamespace export(OdeControl_set_step_size_initial)
NULL
