#' Lorenz Solver R6 Class
#' 
#' @description R6 wrapper for Lorenz ODE solver with optional AD support
#' @export
Lorenz_Solver <- R6::R6Class(
  "Lorenz_Solver",
  public = list(
    ptr = NULL,
    active = FALSE,
    
    initialize = function(System_xp, control_xp, active = FALSE) {
      self$active <- active
      self$ptr <- Solver_new(System_xp, control_xp, active)
    },
    
    time = function() {
      Solver_time(self$ptr, self$active)
    },
    
    state = function() {
      Solver_state(self$ptr, self$active)
    },
    
    set_state = function(y, time) {
      Solver_set_state(self$ptr, y, time, self$active)
      invisible(self)
    },
    
    times = function() {
      Solver_times(self$ptr, self$active)
    },
    
    advance_adaptive = function(times) {
      Solver_advance_adaptive(self$ptr, times, self$active)
      invisible(self)
    },
    
    advance_fixed = function(times) {
      Solver_advance_fixed(self$ptr, times, self$active)
      invisible(self)
    },
    
    step = function() {
      Solver_step(self$ptr, self$active)
      invisible(self)
    },
    
    reset = function() {
      Solver_reset(self$ptr, self$active)
      invisible(self)
    },
    
    collect = function(value) {
      if (missing(value)) {
        Solver_get_collect(self$ptr, self$active)
      } else {
        Solver_set_collect(self$ptr, value, self$active)
        invisible(self)
      }
    },
    
    history_size = function() {
      Solver_get_history_size(self$ptr, self$active)
    },
    
    history_step = function(i) {
      Solver_get_history_step(self$ptr, i, self$active)
    },
    
    history = function() {
      Solver_get_history(self$ptr, self$active) |> 
        dplyr::bind_rows() |> 
        dplyr::as_tibble() |> 
        tibble::remove_rownames()
    },
    
    set_target = function(times, target, obs_indices) {
      Solver_set_target(self$ptr, times, target, obs_indices, self$active)
      invisible(self)
    },
    
    fit = function(ic = NULL, params = NULL) {
      Solver_fit(self$ptr, ic, params)
    }
  )
)

#' Lorenz System R6 Class
#' 
#' @description R6 wrapper for Lorenz system
#' @export
LorenzSystem <- R6::R6Class(
  "LorenzSystem",
  public = list(
    ptr = NULL,
    
    initialize = function(sigma, R, b) {
      self$ptr <- System_new(sigma, R, b)
    },
    
    pars = function() {
      System_pars(self$ptr)
    },
    
    set_params = function(params) {
      System_set_params(self$ptr, params)
      invisible(self)
    }
    
    ,
    set_state = function(y, time = 0.0) {
      System_set_state(self$ptr, y, time)
      invisible(self)
    },
    
    state = function() {
      System_state(self$ptr)
    },
    
    rates = function() {
      System_rates(self$ptr)
    }
    ,
    
    set_initial_state = function(y, t0 = 0.0) {
      System_set_initial_state(self$ptr, y, t0)
      invisible(self)
    },

    reset = function() {
      System_reset(self$ptr)
      invisible(self)
    }
  )
)
