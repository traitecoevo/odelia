#' Lorenz Solver R6 Class
#'
#' @description R6 wrapper for the Lorenz ODE solver.
#' @field ptr External pointer to the underlying C++ solver object.
#' @param System_xp External pointer to the Lorenz system object.
#' @param control_xp External pointer to the ODE control object.
#' @param y Numeric state vector.
#' @param time Scalar time value.
#' @param times Numeric vector of requested output times.
#' @param value Logical flag for collect-history behavior.
#' @param i Integer index into solver history.
#' @param method Integration method: \code{"rkck"} (default, explicit Cash-Karp
#'   RK 4(5)) or \code{"rodas"} (implicit RODAS4(3) Rosenbrock, for stiff
#'   systems).
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

    #' @description Set the calibration target for \code{$fit()}.
    #' @param times Numeric vector: the schedule \code{$fit()} replays, one entry
    #'   per step. Normally a reference run's \code{$times()}, so each fit takes
    #'   the same steps and the loss is a smooth function of the inputs.
    #' @param target Numeric matrix of observed states, one row per entry of
    #'   \code{obs_indices} and one column per state variable.
    #' @param obs_indices Integer indices (1-based) into \code{times} at which
    #'   each row of \code{target} was observed.
    set_target = function(times, target, obs_indices) {
      private$target <- list(times = as.numeric(times),
                             target = as.matrix(target),
                             obs_indices = as.integer(obs_indices))
      invisible(self)
    },

    #' @description Least-squares loss against the target, and its exact
    #'   gradient by reverse-mode automatic differentiation.
    #' @param ic Initial state to fit from, or \code{NULL} to keep the system's.
    #' @param params Parameters to fit at, or \code{NULL} to keep the system's.
    #' @return A list with \code{loss} (sum of squared differences) and
    #'   \code{gradient}:
    #'   d(loss)/d(params), then d(loss)/d(ic), for whichever were given.
    fit = function(ic = NULL, params = NULL) {
      if (is.null(private$target)) {
        stop("Must call set_target() before fit()")
      }
      t <- private$target
      Solver_fit(self$ptr, t$times, t$target, t$obs_indices, ic, params)
    }
  ),
  private = list(
    target = NULL
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
