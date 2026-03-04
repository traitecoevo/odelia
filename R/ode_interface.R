#' OdeControl R6 Class
#' 
#' @description R6 wrapper for ODE solver control settings
#' @export
OdeControl <- R6::R6Class(
  "OdeControl",
  public = list(
    ptr = NULL,
    
    initialize = function() {
      self$ptr <- OdeControl_new()
    },
    
    get_controls = function() {
      OdeControl_get_controls(self$ptr)
    },
    
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
    
    set_tol_abs = function(tol_abs) {
      OdeControl_set_tol_abs(self$ptr, tol_abs)
      invisible(self)
    },
    
    set_tol_rel = function(tol_rel) {
      OdeControl_set_tol_rel(self$ptr, tol_rel)
      invisible(self)
    },
    
    set_a_y = function(a_y) {
      OdeControl_set_a_y(self$ptr, a_y)
      invisible(self)
    },
    
    set_a_dydt = function(a_dydt) {
      OdeControl_set_a_dydt(self$ptr, a_dydt)
      invisible(self)
    },
    
    set_step_size_min = function(step_size_min) {
      OdeControl_set_step_size_min(self$ptr, step_size_min)
      invisible(self)
    },
    
    set_step_size_max = function(step_size_max) {
      OdeControl_set_step_size_max(self$ptr, step_size_max)
      invisible(self)
    },
    
    set_step_size_initial = function(step_size_initial) {
      OdeControl_set_step_size_initial(self$ptr, step_size_initial)
      invisible(self)
    }
  )
)

#' Drivers R6 Class
#' 
#' @description R6 wrapper for external forcing data (drivers)
#' @export
Drivers <- R6::R6Class(
  "Drivers",
  public = list(
    ptr = NULL,
    
    initialize = function() {
      self$ptr <- Drivers_new()
    },
    
    set_constant = function(name, value) {
      Drivers_set_constant(self$ptr, name, value)
      invisible(self)
    },
    
    set_variable = function(name, x, y) {
      Drivers_set_variable(self$ptr, name, x, y)
      invisible(self)
    },
    
    set_extrapolate = function(name, extrapolate) {
      Drivers_set_extrapolate(self$ptr, name, extrapolate)
      invisible(self)
    },
    
    evaluate = function(name, x) {
      Drivers_evaluate(self$ptr, name, x)
    },
    
    evaluate_range = function(name, x) {
      Drivers_evaluate_range(self$ptr, name, x)
    },
    
    get_names = function() {
      Drivers_get_names(self$ptr)
    },
    
    clear = function() {
      Drivers_clear(self$ptr)
      invisible(self)
    }
  )
)