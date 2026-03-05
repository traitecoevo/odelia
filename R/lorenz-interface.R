#' Lorenz Solver R6 Class
#' 
#' @description R6 wrapper for Lorenz ODE solver with optional AD support
#' @field ptr External pointer to the underlying C++ solver object.
#' @field active Logical flag indicating whether AD mode is enabled.
#' @param System_xp External pointer to the Lorenz system object.
#' @param control_xp External pointer to the ODE control object.
#' @param active Logical flag for AD-enabled solver mode.
#' @param y Numeric state vector.
#' @param time Scalar time value.
#' @param times Numeric vector of requested output times.
#' @param value Logical flag for collect-history behavior.
#' @param i Integer index into solver history.
#' @param target Numeric matrix of target observations.
#' @param obs_indices Integer vector of observed-state indices.
#' @param ic Optional initial condition value(s) for fitting.
#' @param params Optional parameter vector for fitting.
#' @export
Lorenz_Solver <- R6::R6Class(
  "Lorenz_Solver",
  public = list(
    ptr = NULL,
    active = FALSE,

    #' @description Initialize a solver for a Lorenz system.
    initialize = function(System_xp, control_xp, active = FALSE) {
      self$active <- active
      self$ptr <- Solver_new(System_xp, control_xp, active)
    },

    #' @description Get current solver time.
    time = function() {
      Solver_time(self$ptr, self$active)
    },

    #' @description Get current solver state.
    state = function() {
      Solver_state(self$ptr, self$active)
    },

    #' @description Set current solver state and time.
    set_state = function(y, time) {
      Solver_set_state(self$ptr, y, time, self$active)
      invisible(self)
    },

    #' @description Get stored solver times.
    times = function() {
      Solver_times(self$ptr, self$active)
    },

    #' @description Advance solver using adaptive stepping.
    advance_adaptive = function(times) {
      Solver_advance_adaptive(self$ptr, times, self$active)
      invisible(self)
    },

    #' @description Advance solver using fixed stepping.
    advance_fixed = function(times) {
      Solver_advance_fixed(self$ptr, times, self$active)
      invisible(self)
    },

    #' @description Advance solver by one step.
    step = function() {
      Solver_step(self$ptr, self$active)
      invisible(self)
    },

    #' @description Reset solver to its initial state.
    reset = function() {
      Solver_reset(self$ptr, self$active)
      invisible(self)
    },

    #' @description Get or set history collection behavior.
    collect = function(value) {
      if (missing(value)) {
        Solver_get_collect(self$ptr, self$active)
      } else {
        Solver_set_collect(self$ptr, value, self$active)
        invisible(self)
      }
    },

    #' @description Return number of stored history entries.
    history_size = function() {
      Solver_get_history_size(self$ptr, self$active)
    },

    #' @description Return one history entry by index.
    history_step = function(i) {
      Solver_get_history_step(self$ptr, i, self$active)
    },

    #' @description Return history as a tibble.
    history = function() {
      Solver_get_history(self$ptr, self$active) |> 
        dplyr::bind_rows() |> 
        dplyr::as_tibble() |> 
        tibble::remove_rownames()
    },

    #' @description Set calibration targets for fitting.
    set_target = function(times, target, obs_indices) {
      Solver_set_target(self$ptr, times, target, obs_indices, self$active)
      invisible(self)
    },

    #' @description Fit initial conditions and/or parameters.
    fit = function(ic = NULL, params = NULL) {
      Solver_fit(self$ptr, ic, params)
    }
  )
)

#' Lorenz System R6 Class
#' 
#' @description R6 wrapper for Lorenz system
#' @field ptr External pointer to the underlying C++ Lorenz system object.
#' @param sigma Lorenz parameter sigma.
#' @param R Lorenz parameter R.
#' @param b Lorenz parameter b.
#' @param params Numeric vector of system parameters.
#' @param y Numeric state vector.
#' @param time Scalar time value.
#' @param t0 Initial time value.
#' @export
LorenzSystem <- R6::R6Class(
  "LorenzSystem",
  public = list(
    ptr = NULL,

    #' @description Initialize a Lorenz system object.
    initialize = function(sigma, R, b) {
      self$ptr <- System_new(sigma, R, b)
    },

    #' @description Return current system parameters.
    pars = function() {
      System_pars(self$ptr)
    },

    #' @description Set model parameters.
    set_params = function(params) {
      System_set_params(self$ptr, params)
      invisible(self)
    }

    ,
    #' @description Set system state and time.
    set_state = function(y, time = 0.0) {
      System_set_state(self$ptr, y, time)
      invisible(self)
    },

    #' @description Return current state.
    state = function() {
      System_state(self$ptr)
    },

    #' @description Return current rates.
    rates = function() {
      System_rates(self$ptr)
    }
    ,

    #' @description Set initial state and initial time.
    set_initial_state = function(y, t0 = 0.0) {
      System_set_initial_state(self$ptr, y, t0)
      invisible(self)
    },

    #' @description Reset the system to its initial condition.
    reset = function() {
      System_reset(self$ptr)
      invisible(self)
    }
  )
)
