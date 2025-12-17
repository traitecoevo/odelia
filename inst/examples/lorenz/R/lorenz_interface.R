Runner <- R6::R6Class(
  "Runner",
  public = list(
    ptr = NULL,
    initialize = function(LorenzSystem_xp, control_xp) {
      self$ptr <- Runner_new(LorenzSystem_xp, control_xp)
    },
    time = function() {
      Runner_time(self$ptr)
    },
    state = function() {
      Runner_state(self$ptr)
    },
    set_state = function(y, time) {
      Runner_set_state(self$ptr, y, time)
      invisible(self)
    },
    times = function() {
      Runner_times(self$ptr)
    },
    advance_adaptive = function(times) {
      Runner_advance_adaptive(self$ptr, times)
      invisible(self)
    },
    advance_fixed = function(times) {
      Runner_advance_fixed(self$ptr, times)
      invisible(self)
    },
    step = function() {
      Runner_step(self$ptr)
      invisible(self)
    },
    reset = function() {
      Runner_reset(self$ptr)
      invisible(self)
    },
    collect = function(value) {
    if (missing(value)) {
      Runner_get_collect(self$ptr)
    } else {
      Runner_set_collect(self$ptr, value)
    }
    },
    history_size = function() {
      Runner_get_history_size(self$ptr)
    },
    history_element = function(i) {
      # 1-based index from R, convert to 0-based size_t if you prefer
      Runner_get_history_element(self$ptr, as.integer(i - 1L))
    },
    history = function() {
      Runner_get_history(self$ptr)
    }
  )
)

LorenzSystem <- R6::R6Class(
  "LorenzSystem",
  public = list(
    ptr = NULL,
    initialize = function(sigma, R, b) {
      # create the C++ LorenzSystem object
      self$ptr <- LorenzSystem_new(sigma, R, b)
    },
    pars = function() {
      LorenzSystem_pars(self$ptr)
    },
    set_state = function(y) {
      LorenzSystem_set_state(self$ptr, y)
      invisible(self)
    },
    state = function() {
      LorenzSystem_state(self$ptr)
    },
    rates = function() {
      LorenzSystem_rates(self$ptr)
    }
  )
)

