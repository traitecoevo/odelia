LeafThermalSolver <- R6::R6Class(
  "LeafThermalSolver",
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
    history_step = function(i) {
      Solver_get_history_step(self$ptr, i)
    },
    history = function() {
      Solver_get_history(self$ptr)|> dplyr::bind_rows() |> dplyr::as_tibble() |> tibble::remove_rownames()
    }
 )
)

LeafThermalSystem <- R6::R6Class(
  "LeafThermalSystem",
  public = list(
    ptr = NULL,
    initialize = function(pars, drivers) {
      # create the C++ LeafThermalSystem object
      self$ptr <- LeafThermalSystem_new(pars, drivers$ptr)
    },
    pars = function() {
      LeafThermalSystem_pars(self$ptr)
    },
    initialize_drivers = function(drivers) {
      LeafThermalSystem_initialize_drivers(self$ptr, drivers$ptr)
      invisible(self)
    },
    set_state = function(y, time) {
      LeafThermalSystem_set_state(self$ptr, y, time)
      invisible(self)
    },
    state = function() {
      LeafThermalSystem_state(self$ptr)
    },
    rates = function() {
      LeafThermalSystem_rates(self$ptr)
    },
    get_current_drivers = function() {
      LeafThermalSystem_get_current_drivers(self$ptr)
    }
  )
)

Drivers <- R6::R6Class(
  "Drivers",
  public = list(
    ptr = NULL,
    initialize = function() {
      # create the C++ Drivers object
      self$ptr <- Drivers_new()
    },
    set_constant = function(driver_name, k) {
      Drivers_set_constant(self$ptr, driver_name, k)
      invisible(self)
    },
    set_variable = function(driver_name, x, y) {
      Drivers_set_variable(self$ptr, driver_name, x, y)
      invisible(self)
    },
    set_extrapolate = function(driver_name, extrapolate) {
      Drivers_set_extrapolate(self$ptr, driver_name, extrapolate)
      invisible(self)
    },
    evaluate = function(driver_name, x) {
      Drivers_evaluate(self$ptr, driver_name, x)
    },
    evaluate_range = function(driver_name, x) {
      Drivers_evaluate_range(self$ptr, driver_name, x)
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
