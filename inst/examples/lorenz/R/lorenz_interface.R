Lorenz_Solver <- R6::R6Class(
  "Lorenz_Solver",
  public = list(
    ptr = NULL,
    initialize = function(System_xp, control_xp) {
      self$ptr <- Solver_new(System_xp, control_xp)
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
    history_element = function(i) {
      # 1-based index from R, convert to 0-based size_t if you prefer
      Solver_get_history_element(self$ptr, as.integer(i - 1L))
    },
    history = function() {
      Solver_get_history(self$ptr)
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
    set_state = function(y) {
      System_set_state(self$ptr, y)
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

