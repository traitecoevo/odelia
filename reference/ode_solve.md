# Solve an ODE given as an R function, deSolve style

A convenience over \[OdeSolver\] shaped like \`deSolve::ode()\`: the
right-hand side is \`func(t, y, parms)\` returning a list whose first
element is the rates (a bare numeric vector is accepted too), the
optional Jacobian is \`jacfunc(t, y, parms)\` returning an n by n matrix
with column j = \`d(dy/dt)/dy\[j\]\`, and the result is a matrix with a
\`time\` column and one column per state variable, one row per requested
time. A function written for deSolve runs unchanged.

## Usage

``` r
ode_solve(
  func,
  y0,
  times,
  parms = NULL,
  jacfunc = NULL,
  method = "dopri",
  control = NULL,
  rtol = 1e-06,
  atol = 1e-06,
  controller = "gsl",
  autonomous = FALSE,
  jac_fd_step = 1e-06,
  jac_fd_floor = 1e-05,
  dense = NULL
)
```

## Arguments

- func:

  The right-hand side, \`func(t, y, parms)\`.

- y0:

  Initial state; names, if any, become column names.

- times:

  Times to report at; the first is the initial time.

- parms:

  Passed through to \`func\` and \`jacfunc\`.

- jacfunc:

  \`NULL\`, or the Jacobian \`jacfunc(t, y, parms)\`.

- method:

  \`"dopri"\` (explicit Dormand–Prince 5(4), the default, as
  \`deSolve\`'s \`ode45\`), \`"rkck"\` (explicit Cash–Karp 4(5)) or
  \`"rodas"\` (implicit RODAS4(3), for stiff problems).

- control:

  An \[OdeControl\], or \`NULL\` for the defaults with the tolerances
  below, the \`controller\` below, and the largest step set to the span
  of \`times\` (an \[OdeControl\] made directly caps the step at 10).
  With a \`control\` given, \`rtol\`, \`atol\` and \`controller\` are
  ignored.

- rtol, atol:

  Shorthand for the control's relative and absolute tolerances when
  \`control\` is \`NULL\`.

- controller:

  The step-size rule when \`control\` is \`NULL\`: \`"gsl"\` (the
  default) or \`"hairer"\`; see \[OdeControl\].

- autonomous:

  \`TRUE\` if \`func\` does not depend on \`t\`.

- jac_fd_step, jac_fd_floor:

  The finite-difference Jacobian's relative step and small-component
  floor; see \[OdeSolver\].

- dense:

  \`TRUE\` to let the stepper choose its own steps and read the
  requested times off the interpolant of the step spanning each, as
  \`deSolve\`'s \`lsoda\` does; \`FALSE\` to land a step on every
  requested time, as its \`ode45\` does, which costs more steps when the
  output is finer than the steps. \`NULL\`, the default, means \`TRUE\`
  under \`"dopri"\`, whose interpolant has the stepper's order, and
  \`FALSE\` under \`"rkck"\` and \`"rodas"\`, whose interpolant is cubic
  Hermite, one order short (on Lorenz at 1e-6 the Cash–Karp dense rows
  err by 6e-4 where landed rows err by 3e-11, measured).

## Value

A numeric matrix of class \`odelia_solution\`, \`length(times)\` rows;
\[ode_counts()\] on it says what the solve cost.

## Examples

``` r
lorenz <- function(t, y, p) {
  list(c(p[["sigma"]] * (y[2] - y[1]),
         y[1] * (p[["rho"]] - y[3]) - y[2],
         y[1] * y[2] - p[["beta"]] * y[3]))
}
out <- ode_solve(lorenz, y0 = c(x = 1, y = 1, z = 1), times = seq(0, 2, by = 0.1),
                 parms = c(sigma = 10, rho = 28, beta = 8 / 3), autonomous = TRUE)
head(out)
#>      time         x         y         z
#> [1,]  0.0  1.000000  1.000000  1.000000
#> [2,]  0.1  2.133108  4.471421  1.113900
#> [3,]  0.2  6.542529 13.731189  4.180196
#> [4,]  0.3 16.684825 27.183509 26.206476
#> [5,]  0.4 15.366195  1.113016 46.757851
#> [6,]  0.5  1.198264 -8.867204 32.454742
ode_counts(out)
#>        n_rhs        n_jac      n_steps n_rejections 
#>          601            0           89           11 
```
