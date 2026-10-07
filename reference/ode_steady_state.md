# Equilibrium of an ODE given as an R function, and how it moves with the parameters

Solve \`func(t, y, parms) = 0\` for a fixed point \`y\*\` by damped
Newton, say whether it is attracting, and take its sensitivity to every
numeric parameter by the implicit function theorem, \`dy\*/dp =
-(df/dy)^-1 (df/dp)\`, with no integration through the transient. The
right-hand side has the shape \[ode_solve()\] takes, so a function
written for \`deSolve\` (or for \`rootSolve::steady()\`) runs unchanged.

## Usage

``` r
ode_steady_state(
  func,
  y0,
  parms = NULL,
  jacfunc = NULL,
  t0 = 0,
  tol = 1e-10,
  max_iter = 100L,
  warmup = NULL,
  method = "rodas",
  control = NULL,
  rtol = 1e-06,
  atol = 1e-06,
  sensitivity = TRUE,
  parms_fd_step = 1e-06,
  autonomous = FALSE,
  jac_fd_step = 1e-06,
  jac_fd_floor = 1e-05,
  line_search = TRUE
)
```

## Arguments

- func:

  The right-hand side, \`func(t, y, parms)\`, returning the rates or a
  list whose first element is the rates.

- y0:

  Starting guess; names, if any, name the result.

- parms:

  Passed through to \`func\` and \`jacfunc\`. Its numeric elements are
  what the sensitivity is taken with respect to.

- jacfunc:

  \`NULL\`, or the Jacobian \`jacfunc(t, y, parms)\`.

- t0:

  The time \`func\` is evaluated at.

- tol:

  Convergence: the largest absolute rate at the solution.

- max_iter:

  The most Newton iterations to take.

- warmup:

  \`NULL\` for Newton alone, or the span of time to integrate the
  transient over before retrying Newton when its root is not attracting;
  a vector of at least two times, the first \`t0\`, is taken as the
  times to land on instead.

- method:

  The stepper for the warm-up: \`"rodas"\` (implicit RODAS4(3), the
  default, for a stiff system), \`"dopri"\` or \`"rkck"\`.

- control:

  An \[OdeControl\] for the warm-up, or \`NULL\` for the defaults with
  \`rtol\` and \`atol\`.

- rtol, atol:

  The warm-up's tolerances when \`control\` is \`NULL\`.

- sensitivity:

  \`TRUE\` for every numeric parameter, \`FALSE\` for none, or the names
  or indices of the parameters to take it for.

- parms_fd_step:

  Relative step of the central difference in each parameter: \`h =
  parms_fd_step \* max(abs(p), 1)\`.

- autonomous:

  \`TRUE\` if \`func\` does not depend on \`t\`, which saves the
  evaluation behind \`time_dependence\`.

- jac_fd_step, jac_fd_floor:

  The finite-difference Jacobian's relative step and small-component
  floor; see \[OdeSolver\].

- line_search:

  \`TRUE\` to backtrack on the residual when a full Newton step would
  not reduce it, which widens the basin of convergence.

## Value

A list of class \`odelia_steady_state\`: \`y\`, the equilibrium estimate
(named as \`y0\`); \`residual\`, the rates there, and \`residual_norm\`,
the largest absolute one; \`iterations\`, \`converged\` and \`warmed\`;
\`jacobian\`, the n by n \`df/dy\` at \`y\*\`, column j the derivative
with respect to \`y\[j\]\`; \`eigenvalues\`, those of the Jacobian,
complex; \`spectral_abscissa\`, the largest real part; and \`stable\`,
\`TRUE\` when that is negative; \`time_dependence\`, the largest
absolute \`df/dt\` at the solution; \`sensitivity\`, the n by p matrix
\`dy\*/dp\`, rows named as \`y\`, columns as the parameters, or \`NULL\`
when none was asked for or possible; and \`counts\`, what the solve cost
(\[ode_counts()\]). Without convergence the Jacobian and everything read
off it are \`NULL\`, and a warning says so.

## Details

Newton needs the Jacobian \`df/dy\`: supply \`jacfunc(t, y, parms)\`
returning the n by n matrix whose column j is the derivative with
respect to \`y\[j\]\`, or leave it \`NULL\` and it is formed by forward
differences at \`n\` right-hand-side evaluations per iteration. Each
iteration also costs one evaluation per line-search trial, usually one.

Newton finds \*a\* root, not necessarily an attracting one: from a guess
near the trivial equilibrium of a demographic model it converges there
in a step or two, and \`stable\` says so. Give \`warmup\` to seek an
attractor instead: when Newton's root is not attracting, or Newton does
not converge, the transient is integrated from \`y0\` for that span with
\`method\` (RODAS for a stiff system) and Newton is retried from where
it ends; \`warmed\` reports that this happened.

The sensitivity is to the parameters \`func\` can be perturbed in from
R: every element of a numeric \`parms\`, or every element of a list
\`parms\` that is a single number. Each column is a central difference
of \`func\` at \`y\*\`, two evaluations per parameter, so its accuracy
is that of the difference (\`parms_fd_step\`), not of the Jacobian.
\`df/dy\` is singular at a bifurcation point, where the theorem does not
apply; the sensitivity stops with an error there.

Equilibrium is only meaningful for an autonomous system.
\`time_dependence\` is the size of \`df/dt\` at the solution, formed by
one extra evaluation unless \`autonomous = TRUE\`; a value above
finite-difference noise means \`func\` depends on \`t\` and the
"equilibrium" is a root at \`t0\` only.

## Examples

``` r
# y0' = a - b y0 ; y1' = y0^2 - c y1 : the fixed point is (a/b, a^2/(b^2 c)).
f <- function(t, y, p) c(p[["a"]] - p[["b"]] * y[1], y[1]^2 - p[["c"]] * y[2])
eq <- ode_steady_state(f, y0 = c(n = 0, m = 0), parms = c(a = 2, b = 1.5, c = 0.7))
eq$y
#>        n        m 
#> 1.333333 2.539683 
eq$stable
#> [1] TRUE
eq$sensitivity   # dy*/da, dy*/db, dy*/dc
#>           a          b         c
#> n 0.6666667 -0.8888889  0.000000
#> m 2.5396838 -3.3862451 -3.628118

# A logistic population has a trivial root at 0 that Newton finds from a
# small guess; warmup integrates past it to the attractor.
g <- function(t, y, p) p[["r"]] * y * (1 - y / p[["K"]])
ode_steady_state(g, y0 = 1e-3, parms = c(r = 1, K = 10))$y
#>           y1 
#> -1.00017e-15 
ode_steady_state(g, y0 = 1e-3, parms = c(r = 1, K = 10), warmup = 50)$y
#> y1 
#> 10 
```
