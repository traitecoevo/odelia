# Lorenz Solver R6 Class

R6 wrapper for the Lorenz ODE solver.

## Public fields

- `ptr`:

  External pointer to the underlying C++ solver object.

## Methods

### Public methods

- [`Lorenz_Solver$new()`](#method-Lorenz_Solver-initialize)

- [`Lorenz_Solver$time()`](#method-Lorenz_Solver-time)

- [`Lorenz_Solver$state()`](#method-Lorenz_Solver-state)

- [`Lorenz_Solver$set_state()`](#method-Lorenz_Solver-set_state)

- [`Lorenz_Solver$times()`](#method-Lorenz_Solver-times)

- [`Lorenz_Solver$advance_adaptive()`](#method-Lorenz_Solver-advance_adaptive)

- [`Lorenz_Solver$advance_fixed()`](#method-Lorenz_Solver-advance_fixed)

- [`Lorenz_Solver$advance_euler()`](#method-Lorenz_Solver-advance_euler)

- [`Lorenz_Solver$step()`](#method-Lorenz_Solver-step)

- [`Lorenz_Solver$reset()`](#method-Lorenz_Solver-reset)

- [`Lorenz_Solver$collect()`](#method-Lorenz_Solver-collect)

- [`Lorenz_Solver$history_size()`](#method-Lorenz_Solver-history_size)

- [`Lorenz_Solver$history_step()`](#method-Lorenz_Solver-history_step)

- [`Lorenz_Solver$history()`](#method-Lorenz_Solver-history)

- [`Lorenz_Solver$set_target()`](#method-Lorenz_Solver-set_target)

- [`Lorenz_Solver$fit()`](#method-Lorenz_Solver-fit)

- [`Lorenz_Solver$clone()`](#method-Lorenz_Solver-clone)

------------------------------------------------------------------------

### `Lorenz_Solver$new()`

Initialize a solver for a Lorenz system.

#### Usage

    Lorenz_Solver$new(System_xp, control_xp, method = "rkck")

#### Arguments

- `System_xp`:

  External pointer to the Lorenz system object.

- `control_xp`:

  External pointer to the ODE control object.

- `method`:

  Integration method: `"rkck"` (default, explicit Cash-Karp RK 4(5)) or
  `"rodas"` (implicit RODAS4(3) Rosenbrock, for stiff systems).

------------------------------------------------------------------------

### `Lorenz_Solver$time()`

Get current solver time.

#### Usage

    Lorenz_Solver$time()

------------------------------------------------------------------------

### `Lorenz_Solver$state()`

Get current solver state.

#### Usage

    Lorenz_Solver$state()

------------------------------------------------------------------------

### `Lorenz_Solver$set_state()`

Set current solver state and time.

#### Usage

    Lorenz_Solver$set_state(y, time)

#### Arguments

- `y`:

  Numeric state vector.

- `time`:

  Scalar time value.

------------------------------------------------------------------------

### `Lorenz_Solver$times()`

Get stored solver times.

#### Usage

    Lorenz_Solver$times()

------------------------------------------------------------------------

### `Lorenz_Solver$advance_adaptive()`

Advance solver using adaptive stepping.

#### Usage

    Lorenz_Solver$advance_adaptive(times)

#### Arguments

- `times`:

  Numeric vector of requested output times.

------------------------------------------------------------------------

### `Lorenz_Solver$advance_fixed()`

Advance solver using fixed stepping.

#### Usage

    Lorenz_Solver$advance_fixed(times)

#### Arguments

- `times`:

  Numeric vector of requested output times.

------------------------------------------------------------------------

### `Lorenz_Solver$advance_euler()`

Advance solver using fixed-step forward Euler.

#### Usage

    Lorenz_Solver$advance_euler(times)

#### Arguments

- `times`:

  Numeric vector of requested output times.

------------------------------------------------------------------------

### `Lorenz_Solver$step()`

Advance solver by one step.

#### Usage

    Lorenz_Solver$step()

------------------------------------------------------------------------

### `Lorenz_Solver$reset()`

Reset solver to its initial state.

#### Usage

    Lorenz_Solver$reset()

------------------------------------------------------------------------

### `Lorenz_Solver$collect()`

Get or set history collection behavior.

#### Usage

    Lorenz_Solver$collect(value)

#### Arguments

- `value`:

  Logical flag for collect-history behavior.

------------------------------------------------------------------------

### `Lorenz_Solver$history_size()`

Return number of stored history entries.

#### Usage

    Lorenz_Solver$history_size()

------------------------------------------------------------------------

### `Lorenz_Solver$history_step()`

Return one history entry by index.

#### Usage

    Lorenz_Solver$history_step(i)

#### Arguments

- `i`:

  Integer index into solver history.

------------------------------------------------------------------------

### `Lorenz_Solver$history()`

Return history as a tibble.

#### Usage

    Lorenz_Solver$history()

------------------------------------------------------------------------

### `Lorenz_Solver$set_target()`

Set the calibration target for `$fit()`.

#### Usage

    Lorenz_Solver$set_target(times, target, obs_indices)

#### Arguments

- `times`:

  Numeric vector: the schedule `$fit()` replays, one entry per step.
  Normally a reference run's `$times()`, so each fit takes the same
  steps and the loss is a smooth function of the inputs.

- `target`:

  Numeric matrix of observed states, one row per entry of `obs_indices`
  and one column per state variable.

- `obs_indices`:

  Integer indices (1-based) into `times` at which each row of `target`
  was observed.

------------------------------------------------------------------------

### `Lorenz_Solver$fit()`

Least-squares loss against the target, and its exact gradient by
reverse-mode automatic differentiation.

#### Usage

    Lorenz_Solver$fit(ic = NULL, params = NULL)

#### Arguments

- `ic`:

  Initial state to fit from, or `NULL` to keep the system's.

- `params`:

  Parameters to fit at, or `NULL` to keep the system's.

#### Returns

A list with `loss` (sum of squared differences) and `gradient`:
d(loss)/d(params), then d(loss)/d(ic), for whichever were given.

------------------------------------------------------------------------

### `Lorenz_Solver$clone()`

The objects of this class are cloneable with this method.

#### Usage

    Lorenz_Solver$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
