#' Lorenz Solver R6 Class
#'
#' @description R6 wrapper for the Lorenz ODE solver. Gradients are computed on
#'   this ordinary (double) solver via `value_and_gradient()`; the AD replay is
#'   built internally, so there is no separate active solver to manage.
#' @field ptr External pointer to the underlying C++ solver object.
#' @param System_xp External pointer to the Lorenz system object.
#' @param control_xp External pointer to the ODE control object.
#' @param y Numeric state vector.
#' @param time Scalar time value.
#' @param times Numeric vector of requested output times.
#' @param value Logical flag for collect-history behavior.
#' @param i Integer index into solver history.
#' @param observations Numeric matrix of measured observations to fit against.
#' @param obs_indices Integer vector of observed-state indices.
#' @param ic Optional initial condition value(s) to differentiate w.r.t.
#' @param params Optional parameter vector to differentiate w.r.t.
#' @param method Integration method: `"rkck"` (default, explicit Cash-Karp
#'   RK 4(5)) or `"rodas"` (implicit RODAS4(3) Rosenbrock, for stiff systems).
#' @export
Lorenz_Solver <- R6::R6Class(
  "Lorenz_Solver",
  public = list(
    ptr = NULL,

    #' @description Initialize a solver for a Lorenz system.
    initialize = function(System_xp, control_xp, method = "rkck") {
      self$ptr <- Solver_new(System_xp, control_xp, method)
    },

    #' @description Get current solver time.
    time = function() {
      Solver_time(self$ptr)
    },

    #' @description Get current solver state.
    state = function() {
      Solver_state(self$ptr)
    },

    #' @description Set current solver state and time.
    set_state = function(y, time) {
      Solver_set_state(self$ptr, y, time)
      invisible(self)
    },

    #' @description Get stored solver times.
    times = function() {
      Solver_times(self$ptr)
    },

    #' @description Advance solver using adaptive stepping.
    advance_adaptive = function(times) {
      Solver_advance_adaptive(self$ptr, times)
      invisible(self)
    },

    #' @description Advance solver using fixed stepping.
    advance_fixed = function(times) {
      Solver_advance_fixed(self$ptr, times)
      invisible(self)
    },

    #' @description Advance solver using fixed-step forward Euler.
    advance_euler = function(times) {
      Solver_advance_euler(self$ptr, times)
      invisible(self)
    },
    #' @description Advance solver by one step.
    step = function() {
      Solver_step(self$ptr)
      invisible(self)
    },

    #' @description Reset solver to its initial state.
    reset = function() {
      Solver_reset(self$ptr)
      invisible(self)
    },

    #' @description Get or set history collection behavior.
    collect = function(value) {
      if (missing(value)) {
        Solver_get_collect(self$ptr)
      } else {
        Solver_set_collect(self$ptr, value)
        invisible(self)
      }
    },

    #' @description Return number of stored history entries.
    history_size = function() {
      Solver_get_history_size(self$ptr)
    },

    #' @description Return one history entry by index.
    history_step = function(i) {
      Solver_get_history_step(self$ptr, i)
    },

    #' @description Return history as a tibble.
    history = function() {
      Solver_get_history(self$ptr) |>
        dplyr::bind_rows() |>
        dplyr::as_tibble() |>
        tibble::remove_rownames()
    },

    #' @description Set the calibration observations to fit against. These are
    #'   held by the R wrapper and passed to each `value_and_gradient()` call;
    #'   the C++ solver stores no calibration state.
    set_observations = function(times, observations, obs_indices) {
      private$obs_times <- times
      private$obs_data <- observations
      private$obs_indices <- obs_indices
      invisible(self)
    },

    #' @description Value and least-squares gradient w.r.t. initial conditions
    #'   and/or parameters, from one AD recording. Requires `set_observations()`.
    value_and_gradient = function(ic = NULL, params = NULL) {
      if (is.null(private$obs_times)) {
        stop("Call set_observations() before value_and_gradient()")
      }
      Solver_value_and_gradient(
        self$ptr, private$obs_times, private$obs_data, private$obs_indices,
        ic, params
      )
    }
  ),
  private = list(
    obs_times = NULL,
    obs_data = NULL,
    obs_indices = NULL
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
