#' OdeControl R6 Class
#' 
#' @description R6 wrapper for ODE solver control settings
#' @field ptr External pointer to the underlying C++ OdeControl object.
#' @param tol_abs Absolute tolerance.
#' @param tol_rel Relative tolerance.
#' @param a_y State scaling coefficient.
#' @param a_dydt Derivative scaling coefficient.
#' @param step_size_min Minimum allowed step size.
#' @param step_size_max Maximum allowed step size.
#' @param step_size_initial Initial step size.
#' @export
OdeControl <- R6::R6Class(
  "OdeControl",
  public = list(
    ptr = NULL,

    #' @description Initialize an `OdeControl` instance.
    initialize = function() {
      self$ptr <- OdeControl_new()
    },

    #' @description Get current solver control settings.
    get_controls = function() {
      OdeControl_get_controls(self$ptr)
    },

    #' @description Set all solver control settings at once.
    set_controls = function(tol_abs, tol_rel,
                            a_y, a_dydt,
                            step_size_min, step_size_max,
                            step_size_initial) {
      OdeControl_set_controls(self$ptr,
                              tol_abs, tol_rel,
                              a_y, a_dydt,
                              step_size_min, step_size_max,
                              step_size_initial)
      invisible(self)
    },

    #' @description Set absolute tolerance.
    set_tol_abs = function(tol_abs) {
      OdeControl_set_tol_abs(self$ptr, tol_abs)
      invisible(self)
    },

    #' @description Set relative tolerance.
    set_tol_rel = function(tol_rel) {
      OdeControl_set_tol_rel(self$ptr, tol_rel)
      invisible(self)
    },

    #' @description Set state scaling coefficient.
    set_a_y = function(a_y) {
      OdeControl_set_a_y(self$ptr, a_y)
      invisible(self)
    },

    #' @description Set derivative scaling coefficient.
    set_a_dydt = function(a_dydt) {
      OdeControl_set_a_dydt(self$ptr, a_dydt)
      invisible(self)
    },

    #' @description Set minimum step size.
    set_step_size_min = function(step_size_min) {
      OdeControl_set_step_size_min(self$ptr, step_size_min)
      invisible(self)
    },

    #' @description Set maximum step size.
    set_step_size_max = function(step_size_max) {
      OdeControl_set_step_size_max(self$ptr, step_size_max)
      invisible(self)
    },

    #' @description Set initial step size.
    set_step_size_initial = function(step_size_initial) {
      OdeControl_set_step_size_initial(self$ptr, step_size_initial)
      invisible(self)
    }
  )
)

#' Drivers R6 Class
#' 
#' @description R6 wrapper for external forcing data (drivers)
#' @field ptr External pointer to the underlying C++ Drivers object.
#' @param name Driver variable name.
#' @param value Constant driver value.
#' @param x Numeric vector of input coordinates (e.g., time).
#' @param y Numeric vector of values corresponding to `x`.
#' @param extrapolate Logical flag controlling extrapolation beyond data range.
#' @export
Drivers <- R6::R6Class(
  "Drivers",
  public = list(
    ptr = NULL,

    #' @description Initialize a `Drivers` instance.
    initialize = function() {
      self$ptr <- Drivers_new()
    },

    #' @description Set a constant driver value.
    set_constant = function(name, value) {
      Drivers_set_constant(self$ptr, name, value)
      invisible(self)
    },

    #' @description Set a variable driver from paired vectors.
    set_variable = function(name, x, y) {
      Drivers_set_variable(self$ptr, name, x, y)
      invisible(self)
    },

    #' @description Configure extrapolation behavior for a driver.
    set_extrapolate = function(name, extrapolate) {
      Drivers_set_extrapolate(self$ptr, name, extrapolate)
      invisible(self)
    },

    #' @description Evaluate a driver at one input value.
    evaluate = function(name, x) {
      Drivers_evaluate(self$ptr, name, x)
    },

    #' @description Evaluate a driver over a vector of input values.
    evaluate_range = function(name, x) {
      Drivers_evaluate_range(self$ptr, name, x)
    },

    #' @description Return names of configured drivers.
    get_names = function() {
      Drivers_get_names(self$ptr)
    },

    #' @description Remove all configured drivers.
    clear = function() {
      Drivers_clear(self$ptr)
      invisible(self)
    }
  )
)
