Lorenz_Solver <- R6::R6Class(
  "Lorenz_Solver",
  public = list(
    ptr = NULL,
    initialize = function(System_xp, control_xp, active = FALSE) {
      self$ptr <- Solver_new(System_xp, control_xp, active)
    },
    time = function() {
      Solver_time(self$ptr)
    },
    state = function() {
      Solver_state(self$ptr)
    },
    set_state = function(y, time) {
      Solver_set_state(self$ptr, y, time)
      invisible(self)
    },
    times = function() {
      Solver_times(self$ptr)
    },
    advance_adaptive = function(times) {
      Solver_advance_adaptive(self$ptr, times)
      invisible(self)
    },
    advance_fixed = function(times) {
      Solver_advance_fixed(self$ptr, times)
      invisible(self)
    },
    step = function() {
      Solver_step(self$ptr)
      invisible(self)
    },
    reset = function() {
      Solver_reset(self$ptr)
      invisible(self)
    },
    collect = function(value) {
    if (missing(value)) {
      Solver_get_collect(self$ptr)
    } else {
      Solver_set_collect(self$ptr, value)
    }
    },
    history_size = function() {
      Solver_get_history_size(self$ptr)
    },
    history_step = function(i) {
      Solver_get_history_step(self$ptr, i)

    },
    history = function() {
      Solver_get_history(self$ptr) |> dplyr::bind_rows() |> dplyr::as_tibble() |> tibble::remove_rownames()
    },
    set_target = function(times, target, obs_indices) {
      Solver_set_target(self$ptr, times, target, obs_indices)
      invisible(self)
    },
    fit_ic = function(ic) {
      Solver_fit_ic(self$ptr, ic)
    },
    fit_params = function(params, ic) {
      Solver_fit_params(self$ptr, params, ic)
    }
  )
)

LorenzSystem <- R6::R6Class(
  "LorenzSystem",
  public = list(
    ptr = NULL,
    initialize = function(sigma, R, b) {
      # create the C++ LorenzSystem object
      self$ptr <- System_new(sigma, R, b)
    },
    pars = function() {
      System_pars(self$ptr)
    },
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
  )
)
#' @export
LorenzSystem

#' @export
Lorenz_Solver
