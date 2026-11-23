# An ODE solver over an R right-hand side

Integrate a system whose right-hand side is an R function, with odelia's
adaptive step control and any of its steppers: the explicit
Dormand–Prince 5(4) pair (the method behind \`deSolve\`'s \`ode45\`,
with dense output of its own order), the explicit Cash–Karp 4(5) pair,
the implicit RODAS4(3) Rosenbrock stepper for stiff problems, or the
ROS34PW2 Rosenbrock-W stepper, which keeps its Jacobian across steps.
The solver is driven a step at a time (or to a sequence of times), and
the state can be reset between steps – at a different length if need be
– which is what a consumer with events needs.

The right-hand side \`rhs(t, y)\` returns the rates as a numeric vector
the length of \`y\`; with \`parms\` given it is called as \`rhs(t, y,
parms)\`, and a list whose first element is the rates is accepted, so a
function written for \`deSolve\` runs unchanged. The implicit steppers
need the Jacobian \`d(dy/dt)/dy\`: supply \`jac(t, y)\` (or \`jac(t, y,
parms)\`) returning the n by n matrix whose column j is the derivative
with respect to \`y\[j\]\`, or leave it \`NULL\` and it is formed by
forward differences at \`n\` extra right-hand-side evaluations per step.

Only the number of right-hand-side evaluations matters for a right-hand
side that is itself expensive, and \`\$counts()\` reports it. Per
accepted step RODAS makes six (five stages and the derivative at the new
point), plus one for the time derivative unless \`autonomous = TRUE\`,
plus the Jacobian once – \`n\` evaluations by finite differences, one
call to \`jac\` otherwise – which is kept across a retry of a rejected
step. A rejected attempt costs its six stage evaluations.

Where the Jacobian is the expensive part, \`"ros34pw2"\` forms it less
often. It is a Rosenbrock-W method of order 3, whose order and error
estimate hold whatever matrix stands in for the Jacobian, so it keeps
one across steps: a Jacobian (and the time derivative with it) serves up
to the control's \`jacobian_max_age\` steps (20 by default; see
\[OdeControl\]), and a new one is formed sooner where a step from a
later point than the one it was formed at is rejected. A step costs four
evaluations, accepted or rejected (three stages and the derivative at
the new point), plus the Jacobian and time derivative on the steps that
form them. It saves most at a loose tolerance: at a tight one it takes
more steps than RODAS, being of order 3, and a kept Jacobian more often
gets a step rejected and is replaced. On Robertson's kinetics to \`t =
1e4\`, with \`atol = 1e-6 \* rtol\` and an analytic Jacobian, it forms 5
Jacobians where RODAS forms 36 at \`rtol = 1e-2\`, 29 against 73 at
\`1e-4\`, and 131 against 187 at \`1e-6\`. A solve under \`"ros34pw2"\`
cannot be differentiated, because the steps at which it formed a
Jacobian are not recorded.

For a cheap right-hand side the cost is instead one R call per stage,
which is why an R right-hand side runs no faster here than through
\`deSolve\`; see the Lorenz benchmark.

A callback that finds its state impossible calls \[domain_error()\] and
the step is rejected and retried smaller. Any other error propagates and
leaves the solver holding a half-finished step; \`set_state()\` makes it
usable again. After every accepted \`step()\`, and after
\`set_state()\`, the last call to \`rhs\` was at exactly \`(time(),
state())\`, so a closure that keeps extra results from its last
evaluation can rely on them matching the solver's state.

## Methods

### Public methods

- [`OdeSolver$new()`](#method-OdeSolver-initialize)

- [`OdeSolver$step()`](#method-OdeSolver-step)

- [`OdeSolver$advance_adaptive()`](#method-OdeSolver-advance_adaptive)

- [`OdeSolver$advance_collect()`](#method-OdeSolver-advance_collect)

- [`OdeSolver$advance_until()`](#method-OdeSolver-advance_until)

- [`OdeSolver$time()`](#method-OdeSolver-time)

- [`OdeSolver$state()`](#method-OdeSolver-state)

- [`OdeSolver$rates()`](#method-OdeSolver-rates)

- [`OdeSolver$times()`](#method-OdeSolver-times)

- [`OdeSolver$set_state()`](#method-OdeSolver-set_state)

- [`OdeSolver$step_size()`](#method-OdeSolver-step_size)

- [`OdeSolver$set_step_size()`](#method-OdeSolver-set_step_size)

- [`OdeSolver$counts()`](#method-OdeSolver-counts)

- [`OdeSolver$clone()`](#method-OdeSolver-clone)

------------------------------------------------------------------------

### `OdeSolver$new()`

Create a solver over an R right-hand side. Evaluates \`rhs\` once, at
\`(t0, y0)\`.

#### Usage

    OdeSolver$new(
      rhs,
      y0,
      t0 = 0,
      jac = NULL,
      state_valid = NULL,
      parms = NULL,
      control = NULL,
      method = "dopri",
      autonomous = FALSE,
      jac_fd_step = 1e-06,
      jac_fd_floor = 1e-05
    )

#### Arguments

- `rhs`:

  Function \`rhs(t, y)\` returning \`dy/dt\`, a numeric vector the
  length of \`y\`.

- `y0`:

  Initial state (numeric, length at least one).

- `t0`:

  Initial time.

- `jac`:

  \`NULL\`, or a function \`jac(t, y)\` returning the n by n Jacobian
  matrix with column j = \`d(dy/dt)/dy\[j\]\`.

- `state_valid`:

  \`NULL\`, or a predicate \`state_valid(t, y)\` returning \`TRUE\` for
  a state the model accepts; a step landing on a refused state is
  rejected.

- `parms`:

  \`NULL\`, or a value passed as a third argument to every callback.

- `control`:

  An \[OdeControl\], or \`NULL\` for the defaults.

- `method`:

  \`"dopri"\` (explicit Dormand–Prince 5(4), the default; the one with
  dense output of its own order), \`"rkck"\` (explicit Cash–Karp 4(5)),
  \`"rodas"\` (implicit RODAS4(3), for stiff problems) or \`"ros34pw2"\`
  (Rosenbrock-W of order 3, for a stiff problem whose Jacobian is
  expensive).

- `autonomous`:

  \`TRUE\` if \`rhs\` does not depend on \`t\`, which saves the implicit
  steppers one evaluation for the time derivative each time they form a
  Jacobian.

- `jac_fd_step`:

  Relative step for the finite-difference Jacobian: \`h_j = jac_fd_step
  \* max(abs(y\[j\]), jac_fd_floor)\`. A right-hand side that is itself
  an iterative solve wants a larger step than the default to stay above
  its own noise.

- `jac_fd_floor`:

  The size below which a state component counts as small for that
  perturbation (the default 1e-5 is \`rodas.f\`'s). A floor of 1
  perturbs a component of size 1e-5 by a tenth of itself; a floor as
  small as a tight absolute tolerance makes the perturbation vanish in
  the subtraction. Lower it only for a state whose components
  legitimately live below 1e-5. A component sitting on the upper edge of
  its domain is perturbed downwards instead when the right-hand side
  refuses the upward point with \[domain_error()\].

------------------------------------------------------------------------

### `OdeSolver$step()`

Take one adaptive step, not passing \`time_max\`. Refused when already
at a finite \`time_max\`.

#### Usage

    OdeSolver$step(time_max = Inf)

#### Arguments

- `time_max`:

  A time the step must not pass; \`Inf\` for no bound.

------------------------------------------------------------------------

### `OdeSolver$advance_adaptive()`

Advance to each of \`times\` in turn by adaptive steps, landing on each
exactly.

#### Usage

    OdeSolver$advance_adaptive(times)

#### Arguments

- `times`:

  Times to advance to; the first must be the current time.

------------------------------------------------------------------------

### `OdeSolver$advance_collect()`

Advance to each of \`times\` in turn and return the state at each: a
matrix with \`time\` in the first column and one column per state
variable, one row per time. The first time must be the current time; the
last is landed on exactly. With \`dense = TRUE\` the steps are the
controller's own and each other time is read off the interpolant of the
step that spans it, so the integration costs the same however many rows
are asked for. Under \`"dopri"\` that interpolant has the stepper's own
order; under the others it is cubic Hermite on the step's endpoints, of
order 3: one order short of \`"rkck"\` and \`"rodas"\`, so prefer
\`"dopri"\` for dense output or \`dense = FALSE\`, which lands a step on
every time, and the order of \`"ros34pw2"\` itself.

#### Usage

    OdeSolver$advance_collect(times, dense = TRUE)

#### Arguments

- `times`:

  Times to advance to; the first must be the current time.

- `dense`:

  \`TRUE\` to read the requested times off each step's interpolant,
  \`FALSE\` to land a step on each; see \`advance_collect()\`.

------------------------------------------------------------------------

### `OdeSolver$advance_until()`

Advance by adaptive steps until \`event(t, y)\` changes the sign it had
when called, or to \`time_max\`, whichever comes first. Returns the time
of the event, or \`NA\` where \`time_max\` came first; either way the
solver is left there. \`event\` returns one number, and a value of zero
at the start is an event at the start.

The crossing is seen at the end of a step and located to \`tol\` in time
on that step's interpolant, at no further evaluation of \`rhs\` (see
\`advance_collect()\` for that interpolant's order under each method).
The solver then steps from the start of that step to the located time,
so the state there is integrated, and it and everything after it are
what a run would do that took the same steps and then landed on that
time. A sign change that starts and ends inside one step is not seen.
Asking again with the same \`event\` starts from a value near zero of
either sign and can stop straight away; ask for the next level instead.

\`event\` is called as \`event(t, y)\`, without \`parms\`; a closure can
carry them. A non-finite value is an error. An error while the crossing
is located leaves the solver at the start of the step that crossed.

#### Usage

    OdeSolver$advance_until(event, time_max, tol = 1e-10)

#### Arguments

- `event`:

  Function \`event(t, y)\` returning one number; see
  \`advance_until()\`.

- `time_max`:

  A time the step must not pass; \`Inf\` for no bound.

- `tol`:

  How closely the time of an event is located: the bracket it is found
  in is narrower than \`tol \* (1 + \|t\|)\`.

------------------------------------------------------------------------

### `OdeSolver$time()`

Current time.

#### Usage

    OdeSolver$time()

------------------------------------------------------------------------

### `OdeSolver$state()`

Current state.

#### Usage

    OdeSolver$state()

------------------------------------------------------------------------

### `OdeSolver$rates()`

Rates at the current state, as last evaluated (no new evaluation after a
step or \`set_state()\`).

#### Usage

    OdeSolver$rates()

------------------------------------------------------------------------

### `OdeSolver$times()`

Times of every accepted step since the last \`set_state()\`, starting
with the time it was set at.

#### Usage

    OdeSolver$times()

------------------------------------------------------------------------

### `OdeSolver$set_state()`

Reset the state and time. \`y\` may have a different length from before.
Resets the step size to the control's initial value and the recorded
times; evaluates \`rhs\` once; clears the effect of an error in a
callback.

#### Usage

    OdeSolver$set_state(y, time)

#### Arguments

- `y`:

  A state vector, of any length.

- `time`:

  The time that state is at.

------------------------------------------------------------------------

### `OdeSolver$step_size()`

The step the controller will try next.

#### Usage

    OdeSolver$step_size()

------------------------------------------------------------------------

### `OdeSolver$set_step_size()`

Set the step the controller tries next, for example to carry a
known-good step across a \`set_state()\`.

#### Usage

    OdeSolver$set_step_size(h)

#### Arguments

- `h`:

  A step size.

------------------------------------------------------------------------

### `OdeSolver$counts()`

What the integration has cost: a named vector with \`n_rhs\`
(right-hand-side evaluations) and \`n_jac\` (Jacobian formations) since
the solver was created, \`n_steps\` (accepted steps since the last
\`set_state()\`) and \`n_rejections\` (attempts rejected and retried
smaller, since the solver was created).

#### Usage

    OdeSolver$counts()

------------------------------------------------------------------------

### `OdeSolver$clone()`

The objects of this class are cloneable with this method.

#### Usage

    OdeSolver$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
decay <- function(t, y) -y
s <- OdeSolver$new(decay, y0 = 1, autonomous = TRUE)
s$advance_adaptive(c(0, 1))
s$state()   # exp(-1)
#> [1] 0.3678794
s$counts()
#>        n_rhs        n_jac      n_steps n_rejections 
#>          109            0           18            0 
```
