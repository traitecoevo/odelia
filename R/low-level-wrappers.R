#' Low-level Drivers and OdeControl wrappers
#'
#' @description
#' Low-level exported wrappers around native C++ interfaces for driver
#' configuration and ODE control. Most users should prefer the higher-level
#' R6 classes [Drivers] and [OdeControl].
#'
#' @name low_level_wrappers
#' @aliases Drivers_new Drivers_set_constant Drivers_set_variable
#' @aliases Drivers_set_extrapolate Drivers_evaluate Drivers_evaluate_range
#' @aliases Drivers_get_names Drivers_clear
#' @aliases OdeControl_new OdeControl_get_controls OdeControl_set_controls
#' @aliases OdeControl_set_tol_abs OdeControl_set_tol_rel OdeControl_set_a_y
#' @aliases OdeControl_set_a_dydt OdeControl_set_step_size_min
#' @aliases OdeControl_set_step_size_max OdeControl_set_step_size_initial
#' @keywords internal
NULL
