Runner <- R6::R6Class(
  "Runner",
  public = list(
    ptr = NULL,

    # LorenzSystem_xp and control_xp are external pointers (XPtr<LorenzSystem>, XPtr<OdeControl>)
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
    set_state_from_system = function() {
      Runner_set_state_from_system(self$ptr)
      invisible(self)
    },
    times = function() {
      Runner_times(self$ptr)
    },
    advance_adaptive = function(time) {
      Runner_advance_adaptive(self$ptr, time)
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
    step_to = function(time) {
      Runner_step_to(self$ptr, time)
      invisible(self)
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

