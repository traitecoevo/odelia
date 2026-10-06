#' Equilibrium of an ODE given as an R function, and how it moves with the parameters
#'
#' Solve `func(t, y, parms) = 0` for a fixed point `y*` by damped Newton, say
#' whether it is attracting, and take its sensitivity to every numeric
#' parameter by the implicit function theorem,
#' `dy*/dp = -(df/dy)^-1 (df/dp)`, with no integration through the transient.
#' The right-hand side has the shape [ode_solve()] takes, so a function
#' written for `deSolve` (or for `rootSolve::steady()`) runs unchanged.
#'
#' Newton needs the Jacobian `df/dy`: supply `jacfunc(t, y, parms)` returning
#' the n by n matrix whose column j is the derivative with respect to `y[j]`,
#' or leave it `NULL` and it is formed by forward differences at `n`
#' right-hand-side evaluations per iteration. Each iteration also costs one
#' evaluation per line-search trial, usually one.
#'
#' Newton finds *a* root, not necessarily an attracting one: from a guess near
#' the trivial equilibrium of a demographic model it converges there in a step
#' or two, and `stable` says so. Give `warmup` to seek an attractor instead:
#' when Newton's root is not attracting, or Newton does not converge, the
#' transient is integrated from `y0` for that span with `method` (RODAS for a
#' stiff system) and Newton is retried from where it ends; `warmed` reports
#' that this happened.
#'
#' The sensitivity is to the parameters `func` can be perturbed in from R: every
#' element of a numeric `parms`, or every element of a list `parms` that is a
#' single number. Each column is a central difference of `func` at `y*`, two
#' evaluations per parameter, so its accuracy is that of the difference
#' (`parms_fd_step`), not of the Jacobian. `df/dy` is singular at a bifurcation
#' point, where the theorem does not apply; the sensitivity stops with an error
#' there.
#'
#' Equilibrium is only meaningful for an autonomous system. `time_dependence`
#' is the size of `df/dt` at the solution, formed by one extra evaluation
#' unless `autonomous = TRUE`; a value above finite-difference noise means
#' `func` depends on `t` and the "equilibrium" is a root at `t0` only.
#'
#' @param func The right-hand side, `func(t, y, parms)`, returning the rates
#'   or a list whose first element is the rates.
#' @param y0 Starting guess; names, if any, name the result.
#' @param parms Passed through to `func` and `jacfunc`. Its numeric elements
#'   are what the sensitivity is taken with respect to.
#' @param jacfunc `NULL`, or the Jacobian `jacfunc(t, y, parms)`.
#' @param t0 The time `func` is evaluated at.
#' @param tol Convergence: the largest absolute rate at the solution.
#' @param max_iter The most Newton iterations to take.
#' @param warmup `NULL` for Newton alone, or the span of time to integrate the
#'   transient over before retrying Newton when its root is not attracting; a
#'   vector of at least two times, the first `t0`, is taken as the times to
#'   land on instead.
#' @param method The stepper for the warm-up: `"rodas"` (implicit RODAS4(3),
#'   the default, for a stiff system), `"dopri"` or `"rkck"`.
#' @param control An [OdeControl] for the warm-up, or `NULL` for the defaults
#'   with `rtol` and `atol`.
#' @param rtol,atol The warm-up's tolerances when `control` is `NULL`.
#' @param sensitivity `TRUE` for every numeric parameter, `FALSE` for none, or
#'   the names or indices of the parameters to take it for.
#' @param parms_fd_step Relative step of the central difference in each
#'   parameter: `h = parms_fd_step * max(abs(p), 1)`.
#' @param autonomous `TRUE` if `func` does not depend on `t`, which saves the
#'   evaluation behind `time_dependence`.
#' @param jac_fd_step,jac_fd_floor The finite-difference Jacobian's relative
#'   step and small-component floor; see [OdeSolver].
#' @param line_search `TRUE` to backtrack on the residual when a full Newton
#'   step would not reduce it, which widens the basin of convergence.
#' @return A list of class `odelia_steady_state`:
#'   `y`, the equilibrium estimate (named as `y0`);
#'   `residual`, the rates there, and `residual_norm`, the largest absolute one;
#'   `iterations`, `converged` and `warmed`;
#'   `jacobian`, the n by n `df/dy` at `y*`, column j the derivative with
#'   respect to `y[j]`;
#'   `eigenvalues`, those of the Jacobian, complex; `spectral_abscissa`, the
#'   largest real part; and `stable`, `TRUE` when that is negative;
#'   `time_dependence`, the largest absolute `df/dt` at the solution;
#'   `sensitivity`, the n by p matrix `dy*/dp`, rows named as `y`, columns
#'   as the parameters, or `NULL` when none was asked for or possible;
#'   and `counts`, what the solve cost ([ode_counts()]).
#'   Without convergence the Jacobian and everything read off it are `NULL`,
#'   and a warning says so.
#' @export
#' @examples
#' # y0' = a - b y0 ; y1' = y0^2 - c y1 : the fixed point is (a/b, a^2/(b^2 c)).
#' f <- function(t, y, p) c(p[["a"]] - p[["b"]] * y[1], y[1]^2 - p[["c"]] * y[2])
#' eq <- ode_steady_state(f, y0 = c(n = 0, m = 0), parms = c(a = 2, b = 1.5, c = 0.7))
#' eq$y
#' eq$stable
#' eq$sensitivity   # dy*/da, dy*/db, dy*/dc
#'
#' # A logistic population has a trivial root at 0 that Newton finds from a
#' # small guess; warmup integrates past it to the attractor.
#' g <- function(t, y, p) p[["r"]] * y * (1 - y / p[["K"]])
#' ode_steady_state(g, y0 = 1e-3, parms = c(r = 1, K = 10))$y
#' ode_steady_state(g, y0 = 1e-3, parms = c(r = 1, K = 10), warmup = 50)$y
ode_steady_state <- function(func, y0, parms = NULL, jacfunc = NULL, t0 = 0,
                             tol = 1e-10, max_iter = 100L, warmup = NULL,
                             method = "rodas", control = NULL,
                             rtol = 1e-6, atol = 1e-6,
                             sensitivity = TRUE, parms_fd_step = 1e-6,
                             autonomous = FALSE, jac_fd_step = 1e-6,
                             jac_fd_floor = 1e-5, line_search = TRUE) {
  func <- check_callback(func, "func")
  jacfunc <- check_callback(jacfunc, "jacfunc", nullable = TRUE)
  # deSolve's form takes three arguments; call it that way even with no parms.
  if (is.null(parms)) parms <- list()
  y_names <- names(y0)
  y0 <- as.numeric(y0)
  if (length(y0) < 1 || anyNA(y0)) {
    stop("y0 must be a numeric vector of at least one finite value", call. = FALSE)
  }
  if (is.null(y_names)) y_names <- paste0("y", seq_along(y0))
  t0 <- as.numeric(t0)
  if (length(t0) != 1 || is.na(t0)) {
    stop("t0 must be a single time", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1 || !(tol > 0)) {
    stop("tol must be a single positive number", call. = FALSE)
  }

  warmup_times <- numeric(0)
  if (!is.null(warmup)) {
    warmup <- as.numeric(warmup)
    if (length(warmup) == 1) {
      if (is.na(warmup) || !(warmup > 0)) {
        stop("warmup must be a positive span of time, or a vector of times", call. = FALSE)
      }
      warmup_times <- c(t0, t0 + warmup)
    } else {
      if (anyNA(warmup) || warmup[1] != t0 || is.unsorted(warmup, strictly = TRUE)) {
        stop("warmup times must start at t0 and increase strictly", call. = FALSE)
      }
      warmup_times <- warmup
    }
  }
  if (is.null(control)) {
    control <- OdeControl$new()
    control$set_tol_rel(rtol)
    control$set_tol_abs(atol)
    # As ode_solve(): the span of the warm-up is the natural bound on a step.
    if (length(warmup_times) >= 2) {
      control$set_step_size_max(diff(range(warmup_times)))
    }
  }
  if (!inherits(control, "OdeControl")) {
    stop("control must be an OdeControl", call. = FALSE)
  }

  res <- RSteadyState_solve(func, jacfunc, parms, y0, t0, isTRUE(autonomous),
                            as.numeric(jac_fd_step), as.numeric(jac_fd_floor),
                            tol, as.integer(max_iter), isTRUE(line_search),
                            1e-10, warmup_times, control$ptr, method)

  names(res$y) <- y_names
  names(res$residual) <- y_names
  n_rhs <- res$n_rhs
  n_jac <- res$n_jac
  res$n_rhs <- NULL
  res$n_jac <- NULL
  if (!is.null(res$jacobian)) {
    dimnames(res$jacobian) <- list(y_names, y_names)
  }

  res$sensitivity <- NULL
  if (!res$converged) {
    warning(sprintf("Newton did not converge in %d iterations (largest rate %.3g)",
                    res$iterations, res$residual_norm), call. = FALSE)
  } else {
    idx <- steady_state_parameters(parms, sensitivity)
    if (length(idx) > 0) {
      res$sensitivity <- steady_state_sensitivity(func, t0, res$y, parms, idx,
                                                  res$jacobian, parms_fd_step)
      n_rhs <- n_rhs + 2 * length(idx)
    }
  }
  res$counts <- c(n_rhs = n_rhs, n_jac = n_jac, n_iterations = res$iterations)
  class(res) <- "odelia_steady_state"
  res
}

# Which elements of parms the sensitivity is taken for: positions into parms,
# named where parms is. `which` is TRUE (every numeric element), FALSE (none),
# or names or indices to keep.
steady_state_parameters <- function(parms, which) {
  if (isFALSE(which)) {
    return(integer(0))
  }
  if (is.numeric(parms)) {
    idx <- seq_along(parms)
  } else if (is.list(parms)) {
    idx <- which(vapply(parms, function(p) is.numeric(p) && length(p) == 1, logical(1)))
  } else {
    idx <- integer(0)
  }
  nm <- names(parms)
  if (!is.null(nm)) names(idx) <- nm[idx]
  if (isTRUE(which)) {
    return(idx)
  }
  if (is.character(which)) {
    missing <- setdiff(which, names(idx))
    if (length(missing) > 0) {
      stop("sensitivity names not among the numeric parameters: ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    return(idx[which])
  }
  if (is.numeric(which)) {
    which <- as.integer(which)
    if (!all(which %in% idx)) {
      stop("sensitivity indices must point at numeric elements of parms", call. = FALSE)
    }
    return(idx[match(which, idx)])
  }
  stop("sensitivity must be TRUE, FALSE, or names or indices into parms", call. = FALSE)
}

# dy*/dp by the implicit function theorem, with df/dp a central difference of
# func in each chosen parameter at y*.
steady_state_sensitivity <- function(func, t0, y, parms, idx, J, rel_step) {
  n <- length(y)
  rates <- function(p) {
    out <- func(t0, y, p)
    if (is.list(out)) out <- out[[1]]
    out <- as.numeric(out)
    if (length(out) != n) {
      stop(sprintf("func returned %d rates, expected %d", length(out), n), call. = FALSE)
    }
    out
  }
  dfdp <- matrix(0, n, length(idx))
  for (k in seq_along(idx)) {
    j <- idx[[k]]
    p0 <- parms[[j]]
    h <- rel_step * max(abs(p0), 1)
    up <- parms
    up[[j]] <- p0 + h
    down <- parms
    down[[j]] <- p0 - h
    dfdp[, k] <- (rates(up) - rates(down)) / (2 * h)
  }
  S <- tryCatch(solve(J, -dfdp), error = function(e) {
    stop("df/dy is singular at the equilibrium (a bifurcation point): ",
         "the implicit function theorem does not give a sensitivity (",
         conditionMessage(e), ")", call. = FALSE)
  })
  cn <- names(idx)
  if (is.null(cn) || any(is.na(cn) | cn == "")) cn <- paste0("p", unname(idx))
  dimnames(S) <- list(names(y), cn)
  S
}

#' @rdname ode_counts
#' @export
ode_counts.odelia_steady_state <- function(x, ...) x$counts

#' @export
print.odelia_steady_state <- function(x, ...) {
  cat("<odelia steady state:",
      if (x$converged) "converged" else "NOT converged",
      sprintf("in %d iterations", x$iterations),
      if (isTRUE(x$warmed)) "after a warm-up" else "",
      ">\n")
  print(x$y, ...)
  if (x$converged) {
    cat(sprintf("stable: %s (spectral abscissa %.3g)\n",
                if (isTRUE(x$stable)) "yes" else "no", x$spectral_abscissa))
    if (!is.null(x$sensitivity)) {
      cat("sensitivity dy*/dp:\n")
      print(x$sensitivity, ...)
    }
  }
  invisible(x)
}
