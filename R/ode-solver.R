#' Signal that a state is outside the model's domain
#'
#' Call this from a right-hand side, Jacobian or validity predicate handed to
#' [OdeSolver] or [ode_solve()], at any depth below it, to say that the state
#' it was given is one the model has no meaning for -- a pool below zero, a
#' probability above one, an inner solve that did not converge -- as distinct
#' from an error in the code. The callback is left at once, as by `return()`
#' (its `on.exit()` expressions run), and the solver *rejects the current
#' step* and retries it smaller, exactly as it does for a compiled system
#' throwing `util::DomainError`. Any other error raised in a callback
#' propagates unchanged and ends the solve, which is what a bug should do.
#'
#' If the smallest permitted step still leaves the domain, the solve stops with
#' an error naming this message, so make it say what was wrong. Called outside
#' a solver callback, this is an ordinary error of class `odelia_domain_error`.
#'
#' @param message What is wrong with the state.
#' @return Does not return.
#' @export
#' @examples
#' rhs <- function(t, y) {
#'   if (y[1] < 0) domain_error("y must stay non-negative")
#'   -sqrt(y[1])
#' }
domain_error <- function(message) {
  message <- as.character(message)
  # The adapter marks its call into the callback; leave that frame with the
  # sentinel the adapter turns into a step rejection. See src/r_system.h.
  calls <- sys.calls()
  for (k in rev(seq_along(calls))) {
    if (isTRUE(attr(calls[[k]], "odelia_callback"))) {
      odelia_return_from(sys.frame(k),
                         structure(NA_real_, odelia_domain_error = message))
    }
  }
  cond <- structure(
    class = c("odelia_domain_error", "error", "condition"),
    list(message = message, call = sys.call(-1))
  )
  stop(cond)
}

check_callback <- function(f, name, nullable = FALSE) {
  if (is.null(f) && nullable) {
    return(NULL)
  }
  if (!is.function(f)) {
    stop(sprintf("%s must be a function of (t, y)", name), call. = FALSE)
  }
  f
}

#' An ODE solver over an R right-hand side
#'
#' @description
#' Integrate a system whose right-hand side is an R function, with odelia's
#' adaptive step control and any of its steppers: the explicit Dormand--Prince
#' 5(4) pair (the method behind `deSolve`'s `ode45`, with dense output of its
#' own order), the explicit Cash--Karp 4(5) pair, or the implicit RODAS4(3)
#' Rosenbrock stepper for stiff problems. The solver is
#' driven a step at a time (or to a sequence of times), and the state can be
#' re-seeded between steps -- at a different length if need be -- which is
#' what a consumer with events needs.
#'
#' The right-hand side `rhs(t, y)` returns the rates as a numeric vector the
#' length of `y`; with `parms` given it is called as `rhs(t, y, parms)`, and
#' a list whose first element is the rates is accepted, so a function written
#' for `deSolve` runs unchanged. The implicit stepper needs the Jacobian
#' `d(dy/dt)/dy`: supply `jac(t, y)` (or `jac(t, y, parms)`) returning the n
#' by n matrix whose column j is the derivative with respect to `y[j]`, or
#' leave it `NULL` and it is formed by forward differences at `n` extra
#' right-hand-side evaluations per step.
#'
#' Only the number of right-hand-side evaluations matters for a right-hand
#' side that is itself expensive, and `counts()` reports it. Per accepted
#' step RODAS makes six (five stages and the derivative at the new point),
#' plus one for the time derivative unless `autonomous = TRUE`, plus the
#' Jacobian once -- `n` evaluations by finite differences, one call to `jac`
#' otherwise -- which is kept across a retry of a rejected step. A rejected
#' attempt costs its six stage evaluations. For a cheap right-hand side the
#' cost is instead one R call per stage, which is why an R right-hand side
#' runs no faster here than through `deSolve`; see the Lorenz benchmark.
#'
#' A callback that finds its state impossible calls [domain_error()] and the
#' step is rejected and retried smaller. Any other error propagates and leaves
#' the solver holding a half-finished step; `set_state()` makes it usable
#' again. After every accepted `step()`, and after `set_state()`, the last
#' call to `rhs` was at exactly `(time(), state())`, so a closure that keeps
#' extra results from its last evaluation can rely on them matching the
#' solver's state.
#'
#' @param rhs Function `rhs(t, y)` returning `dy/dt`, a numeric vector the
#'   length of `y`.
#' @param y0 Initial state (numeric, length at least one).
#' @param t0 Initial time.
#' @param jac `NULL`, or a function `jac(t, y)` returning the n by n Jacobian
#'   matrix with column j = `d(dy/dt)/dy[j]`.
#' @param state_valid `NULL`, or a predicate `state_valid(t, y)` returning
#'   `TRUE` for a state the model accepts; a step landing on a refused state is
#'   rejected.
#' @param parms `NULL`, or a value passed as a third argument to every
#'   callback.
#' @param control An [OdeControl], or `NULL` for the defaults.
#' @param method `"dopri"` (explicit Dormand--Prince 5(4), the default; the
#'   one with dense output of its own order), `"rkck"` (explicit Cash--Karp
#'   4(5)) or `"rodas"` (implicit RODAS4(3), for stiff problems).
#' @param autonomous `TRUE` if `rhs` does not depend on `t`, which saves the
#'   implicit stepper one evaluation per step for the time derivative.
#' @param jac_fd_step Relative step for the finite-difference Jacobian:
#'   `h_j = jac_fd_step * max(abs(y[j]), 1)`. A right-hand side that is itself
#'   an iterative solve wants a larger step than the default to stay above its
#'   own noise.
#' @param time_max A time the step must not pass; `Inf` for no bound.
#' @param times Times to advance to; the first must be the current time.
#' @param y A state vector, of any length.
#' @param time The time that state is at.
#' @param h A step size.
#' @export
#' @examples
#' decay <- function(t, y) -y
#' s <- OdeSolver$new(decay, y0 = 1, autonomous = TRUE)
#' s$advance_adaptive(c(0, 1))
#' s$state()   # exp(-1)
#' s$counts()
OdeSolver <- R6::R6Class(
  "OdeSolver",
  public = list(
    #' @description Create a solver over an R right-hand side. Evaluates
    #'   `rhs` once, at `(t0, y0)`.
    initialize = function(rhs, y0, t0 = 0, jac = NULL, state_valid = NULL,
                          parms = NULL, control = NULL, method = "dopri",
                          autonomous = FALSE, jac_fd_step = 1e-6) {
      rhs <- check_callback(rhs, "rhs")
      jac <- check_callback(jac, "jac", nullable = TRUE)
      state_valid <- check_callback(state_valid, "state_valid", nullable = TRUE)
      y0 <- as.numeric(y0)
      if (length(y0) < 1 || anyNA(y0)) {
        stop("y0 must be a numeric vector of at least one finite value", call. = FALSE)
      }
      if (is.null(control)) {
        control <- OdeControl$new()
      }
      if (!inherits(control, "OdeControl")) {
        stop("control must be an OdeControl", call. = FALSE)
      }
      private$control <- control
      private$ptr <- RSolver_new(rhs, jac, state_valid, parms, y0,
                                 as.numeric(t0), control$ptr, method,
                                 isTRUE(autonomous), jac_fd_step)
    },

    #' @description Take one adaptive step, not passing `time_max`. Refused
    #'   when already at a finite `time_max`.
    step = function(time_max = Inf) {
      private$guard()
      RSolver_step(private$ptr, time_max)
      invisible(self)
    },

    #' @description Advance to each of `times` in turn by adaptive steps,
    #'   landing on each exactly.
    advance_adaptive = function(times) {
      private$guard()
      RSolver_advance_adaptive(private$ptr, as.numeric(times))
      invisible(self)
    },

    #' @description Advance to each of `times` in turn and return the state
    #'   at each: a matrix with `time` in the first column and one column per
    #'   state variable, one row per time. The first time must be the current
    #'   time; the last is landed on exactly. With `dense = TRUE` the steps
    #'   are the controller's own and each other time is read off the
    #'   interpolant of the step that spans it, so the integration costs the
    #'   same however many rows are asked for. Under `"dopri"` that
    #'   interpolant has the stepper's own order; under the other two it is
    #'   cubic Hermite on the step's endpoints, one order short, so prefer
    #'   `"dopri"` for dense output or `dense = FALSE`, which lands a step on
    #'   every time.
    advance_collect = function(times, dense = TRUE) {
      private$guard()
      RSolver_advance_collect(private$ptr, as.numeric(times), isTRUE(dense))
    },

    #' @description Current time.
    time = function() RSolver_time(private$ptr),

    #' @description Current state.
    state = function() RSolver_state(private$ptr),

    #' @description Rates at the current state, as last evaluated (no new
    #'   evaluation after a step or `set_state()`).
    rates = function() RSolver_rates(private$ptr),

    #' @description Times of every accepted step since the last
    #'   `set_state()`, starting with the time it was seeded at.
    times = function() RSolver_times(private$ptr),

    #' @description Re-seed the state and time. `y` may have a different
    #'   length from before. Resets the step size to the control's initial
    #'   value and the recorded times; evaluates `rhs` once; clears the
    #'   effect of an error in a callback.
    set_state = function(y, time) {
      y <- as.numeric(y)
      if (length(y) < 1 || anyNA(y)) {
        stop("y must be a numeric vector of at least one finite value", call. = FALSE)
      }
      RSolver_set_state(private$ptr, y, as.numeric(time))
      invisible(self)
    },

    #' @description The step the controller will try next.
    step_size = function() RSolver_step_size(private$ptr),

    #' @description Set the step the controller tries next, for example to
    #'   carry a known-good step across a `set_state()`.
    set_step_size = function(h) {
      RSolver_set_step_size(private$ptr, as.numeric(h))
      invisible(self)
    },

    #' @description What the integration has cost: a named vector with
    #'   `n_rhs` (right-hand-side evaluations) and `n_jac` (Jacobian
    #'   formations) since the solver was created, `n_steps` (accepted steps
    #'   since the last `set_state()`) and `n_rejections` (attempts rejected
    #'   and retried smaller, since the solver was created).
    counts = function() RSolver_counts(private$ptr)
  ),
  private = list(
    ptr = NULL,
    control = NULL,
    # The solver itself knows whether an attempt was abandoned by an error
    # that escaped from a callback; an error raised before any step (bad
    # arguments) leaves it usable.
    guard = function() {
      if (RSolver_mid_step(private$ptr)) {
        stop("the previous step was interrupted by an error in a callback; ",
             "call set_state() before stepping again", call. = FALSE)
      }
    }
  )
)

#' Solve an ODE given as an R function, deSolve style
#'
#' A convenience over [OdeSolver] shaped like `deSolve::ode()`: the right-hand
#' side is `func(t, y, parms)` returning a list whose first element is the
#' rates (a bare numeric vector is accepted too), the optional Jacobian is
#' `jacfunc(t, y, parms)` returning an n by n matrix with column j =
#' `d(dy/dt)/dy[j]`, and the result is a matrix with a `time` column and one
#' column per state variable, one row per requested time. A function written
#' for deSolve runs unchanged.
#'
#' @param func The right-hand side, `func(t, y, parms)`.
#' @param y0 Initial state; names, if any, become column names.
#' @param times Times to report at; the first is the initial time.
#' @param parms Passed through to `func` and `jacfunc`.
#' @param jacfunc `NULL`, or the Jacobian `jacfunc(t, y, parms)`.
#' @param method `"dopri"` (explicit Dormand--Prince 5(4), the default, as
#'   `deSolve`'s `ode45`), `"rkck"` (explicit Cash--Karp 4(5)) or `"rodas"`
#'   (implicit RODAS4(3), for stiff problems).
#' @param control An [OdeControl], or `NULL` for the defaults.
#' @param rtol,atol Shorthand for the control's relative and absolute
#'   tolerances when `control` is `NULL`.
#' @param autonomous `TRUE` if `func` does not depend on `t`.
#' @param jac_fd_step Relative step for the finite-difference Jacobian; see
#'   [OdeSolver].
#' @param dense `TRUE` (the default) to let the stepper choose its own steps
#'   and read the requested times off the interpolant of the step spanning
#'   each, as `deSolve`'s `lsoda` does; `FALSE` to land a step on every
#'   requested time, as its `ode45` does, which costs more steps when the
#'   output is finer than the steps. Under `"dopri"` the interpolant has the
#'   stepper's order; under `"rkck"` and `"rodas"` it is cubic Hermite, one
#'   order short, so use `dense = FALSE` with those when the output grid is
#'   coarser than the steps.
#' @return A numeric matrix of class `odelia_solution`, `length(times)` rows;
#'   [counts()] on it says what the solve cost.
#' @export
#' @examples
#' lorenz <- function(t, y, p) {
#'   list(c(p[["sigma"]] * (y[2] - y[1]),
#'          y[1] * (p[["rho"]] - y[3]) - y[2],
#'          y[1] * y[2] - p[["beta"]] * y[3]))
#' }
#' out <- ode_solve(lorenz, y0 = c(x = 1, y = 1, z = 1), times = seq(0, 2, by = 0.1),
#'                  parms = c(sigma = 10, rho = 28, beta = 8 / 3), autonomous = TRUE)
#' head(out)
#' counts(out)
ode_solve <- function(func, y0, times, parms = NULL, jacfunc = NULL,
                      method = "dopri", control = NULL,
                      rtol = 1e-6, atol = 1e-6,
                      autonomous = FALSE, jac_fd_step = 1e-6, dense = TRUE) {
  func <- check_callback(func, "func")
  # deSolve's form takes three arguments; call it that way even with no parms.
  if (is.null(parms)) parms <- list()
  times <- as.numeric(times)
  if (length(times) < 1 || anyNA(times) || is.unsorted(times, strictly = TRUE)) {
    stop("times must be strictly increasing", call. = FALSE)
  }
  if (is.null(control)) {
    control <- OdeControl$new()
    control$set_tol_rel(rtol)
    control$set_tol_abs(atol)
  }
  # func and jacfunc are called as f(t, y, parms) and a list result is
  # unwrapped, both inside the adapter: no R closure sits between the solver
  # and the user's function.
  s <- OdeSolver$new(func, y0, t0 = times[1], jac = jacfunc, parms = parms,
                     control = control, method = method,
                     autonomous = autonomous, jac_fd_step = jac_fd_step)
  out <- s$advance_collect(times, dense = dense)
  nm <- names(y0)
  if (is.null(nm)) nm <- paste0("y", seq_along(y0))
  colnames(out) <- c("time", nm)
  attr(out, "counts") <- s$counts()
  class(out) <- c("odelia_solution", class(out))
  out
}

#' What a solve cost
#'
#' The number of right-hand-side evaluations (`n_rhs`), Jacobian formations
#' (`n_jac`), accepted steps (`n_steps`) and rejected attempts
#' (`n_rejections`) behind a result of [ode_solve()], as a named vector. For
#' a right-hand side that is itself expensive these are the whole cost of the
#' integration; `n_rejections` is what the step-size controller wasted.
#'
#' @param x A result of [ode_solve()].
#' @param ... Ignored.
#' @return A named numeric vector.
#' @export
#' @examples
#' out <- ode_solve(function(t, y, p) -y, y0 = 1, times = c(0, 1))
#' counts(out)
counts <- function(x, ...) UseMethod("counts")

#' @rdname counts
#' @export
counts.odelia_solution <- function(x, ...) attr(x, "counts")

#' @export
print.odelia_solution <- function(x, ...) {
  print(unclass(x)[, , drop = FALSE], ...)
  n <- attr(x, "counts")
  cat(sprintf("<%d steps, %d evaluations; counts() for detail>\n",
              as.integer(n[["n_steps"]]), as.integer(n[["n_rhs"]])))
  invisible(x)
}
