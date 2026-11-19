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
#' @param params Numeric vector of system parameters, one entry per parameter
#'   in the system's order (for Lorenz: sigma, R, b).
#' @param method Integration method: \code{"rkck"} (default, explicit Cash-Karp
#'   RK 4(5)), \code{"dopri"} (explicit Dormand-Prince 5(4)) or \code{"rodas"}
#'   (implicit RODAS4(3) Rosenbrock, for stiff systems).
#' @export
Lorenz_Solver <- R6::R6Class(
  "Lorenz_Solver",
  public = list(
    ptr = NULL,

    #' @description Initialize a solver for a Lorenz system. The solver takes
    #'   its own copy of the system as it is at this moment and steps that copy:
    #'   \code{set_params()} or \code{set_state()} on the \code{LorenzSystem}
    #'   afterwards do not reach the solver. Use the solver's own
    #'   \code{$set_params()} and \code{$set_state()}.
    initialize = function(System_xp, control_xp, method = "rkck") {
      self$ptr <- Solver_new(System_xp, control_xp, method)
    },

    #' @description Parameters of the system this solver steps.
    pars = function() {
      Solver_params(self$ptr)
    },

    #' @description Set the parameters this solver steps with, and re-seed it
    #'   at its current state and time as \code{$set_state()} would. Call
    #'   \code{$reset()} as well to start again from the initial state.
    set_params = function(params) {
      Solver_set_params(self$ptr, params)
      invisible(self)
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

    #' @description Set the calibration target for
    #'   \code{$loss_and_gradient()}.
    #' @param schedule Numeric vector: the step times
    #'   \code{$loss_and_gradient()} replays, one entry per step. Normally a
    #'   reference run's \code{$times()}, so each evaluation takes the same steps
    #'   and the loss is a smooth function of the inputs.
    #' @param target Numeric matrix of observed states, one row per entry of
    #'   \code{obs_indices} and one column per state variable.
    #' @param obs_indices Integer indices (1-based) into \code{schedule} at which
    #'   each row of \code{target} was observed.
    set_target = function(schedule, target, obs_indices) {
      private$target <- list(schedule = as.numeric(schedule),
                             target = as.matrix(target),
                             obs_indices = as.integer(obs_indices))
      invisible(self)
    },

    #' @description Least-squares loss against the target, and its exact
    #'   gradient by automatic differentiation. The schedule is stepped under
    #'   this solver's \code{method}.
    #' @param initial_state Initial state to evaluate at, or \code{NULL} to keep
    #'   the system's.
    #' @param params Parameters to evaluate at, or \code{NULL} to keep the
    #'   system's.
    #' @param mode \code{"reverse"}: one reverse sweep, whose cost does not grow
    #'   with the number of inputs. \code{"forward"}: the schedule walked again
    #'   carrying one direction per input, whose cost does; it is the only mode
    #'   under \code{"rodas"}, which cannot be swept. \code{NULL} (the default)
    #'   is \code{"forward"} for a \code{"rodas"} solver and \code{"reverse"}
    #'   otherwise.
    #' @return A list with \code{loss} (sum of squared differences) and
    #'   \code{gradient}: d(loss)/d(params), then d(loss)/d(initial_state), for
    #'   whichever were given.
    loss_and_gradient = function(initial_state = NULL, params = NULL,
                                 mode = NULL) {
      if (is.null(private$target)) {
        stop("Must call set_target() before loss_and_gradient()")
      }
      mode <- if (is.null(mode)) "" else match.arg(mode, c("reverse", "forward"))
      t <- private$target
      Solver_loss_and_gradient(self$ptr, t$schedule, t$target, t$obs_indices,
                               initial_state, params, mode)
    },

    #' @description The state at each of \code{times} and its derivative with
    #'   respect to the parameters and to the state at \code{times[1]}, by
    #'   forward-mode automatic differentiation. The solve is the one
    #'   \code{$advance_adaptive(times)} would take from where this solver
    #'   stands, with a step landing on every time; the solver itself does not
    #'   move. Those steps are then walked again carrying one direction per
    #'   input, under the solver's own \code{method}, \code{"rodas"} included.
    #'   The cost grows with the number of inputs and not with the number of
    #'   times.
    #' @param parameters The parameters to differentiate with respect to:
    #'   \code{TRUE} (all), \code{FALSE} (none), or names (\code{"sigma"},
    #'   \code{"R"}, \code{"b"}) or indices.
    #' @param initial The components of the state at \code{times[1]} to
    #'   differentiate with respect to: \code{TRUE}, \code{FALSE}, or names
    #'   (\code{"x"}, \code{"y"}, \code{"z"}) or indices.
    #' @return A list with \code{time}; \code{y}, the state, one row per time;
    #'   and \code{dy}, an array indexed by time, state component and input:
    #'   the parameters first, then the initial state, named \code{x0},
    #'   \code{y0}, \code{z0}.
    sensitivities = function(times, parameters = TRUE, initial = TRUE) {
      states <- c("x", "y", "z")
      p <- forward_directions(parameters, c("sigma", "R", "b"), "parameters")
      s <- forward_directions(initial, states, "initial")
      out <- Solver_sensitivities(self$ptr, as.numeric(times), p - 1L, s - 1L)
      colnames(out$y) <- states
      dimnames(out$dy) <- list(NULL, states,
                               c(names(p), sprintf("%s0", names(s))))
      out
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

# The inputs a forward pass carries a direction for: positions into `names`
# (1-based), named. `which` is TRUE (all), FALSE (none), or names or indices.
forward_directions <- function(which, names, what) {
  idx <- seq_along(names)
  names(idx) <- names
  if (isTRUE(which)) {
    return(idx)
  }
  if (isFALSE(which) || length(which) == 0) {
    return(idx[0])
  }
  if (is.character(which)) {
    missing <- setdiff(which, names)
    if (length(missing) > 0) {
      stop("'", what, "' names not among ", paste(names, collapse = ", "), ": ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    return(idx[which])
  }
  if (is.numeric(which)) {
    which <- as.integer(which)
    if (anyNA(which) || any(which < 1L | which > length(names))) {
      stop("'", what, "' indices must lie in 1..", length(names), call. = FALSE)
    }
    return(idx[which])
  }
  stop("'", what, "' must be TRUE, FALSE, or names or indices", call. = FALSE)
}
