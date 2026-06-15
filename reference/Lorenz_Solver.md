# Lorenz Solver R6 Class

R6 wrapper for Lorenz ODE solver with optional AD support

## Public fields

- `ptr`:

  External pointer to the underlying C++ solver object.

- `active`:

  Logical flag indicating whether AD mode is enabled.

## Methods

### Public methods

- [`Lorenz_Solver$new()`](#method-Lorenz_Solver-new)

- [`Lorenz_Solver$time()`](#method-Lorenz_Solver-time)

- [`Lorenz_Solver$state()`](#method-Lorenz_Solver-state)

- [`Lorenz_Solver$set_state()`](#method-Lorenz_Solver-set_state)

- [`Lorenz_Solver$times()`](#method-Lorenz_Solver-times)

- [`Lorenz_Solver$advance_adaptive()`](#method-Lorenz_Solver-advance_adaptive)

- [`Lorenz_Solver$advance_fixed()`](#method-Lorenz_Solver-advance_fixed)

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

### Method `new()`

Initialize a solver for a Lorenz system.

#### Usage

    Lorenz_Solver$new(System_xp, control_xp, active = FALSE)

#### Arguments

- `System_xp`:

  External pointer to the Lorenz system object.

- `control_xp`:

  External pointer to the ODE control object.

- `active`:

  Logical flag for AD-enabled solver mode.

------------------------------------------------------------------------

### Method [`time()`](https://rdrr.io/r/stats/time.html)

Get current solver time.

#### Usage

    Lorenz_Solver$time()

------------------------------------------------------------------------

### Method `state()`

Get current solver state.

#### Usage

    Lorenz_Solver$state()

------------------------------------------------------------------------

### Method `set_state()`

Set current solver state and time.

#### Usage

    Lorenz_Solver$set_state(y, time)

#### Arguments

- `y`:

  Numeric state vector.

- `time`:

  Scalar time value.

------------------------------------------------------------------------

### Method `times()`

Get stored solver times.

#### Usage

    Lorenz_Solver$times()

------------------------------------------------------------------------

### Method `advance_adaptive()`

Advance solver using adaptive stepping.

#### Usage

    Lorenz_Solver$advance_adaptive(times)

#### Arguments

- `times`:

  Numeric vector of requested output times.

------------------------------------------------------------------------

### Method `advance_fixed()`

Advance solver using fixed stepping.

#### Usage

    Lorenz_Solver$advance_fixed(times)

#### Arguments

- `times`:

  Numeric vector of requested output times.

------------------------------------------------------------------------

### Method [`step()`](https://rdrr.io/r/stats/step.html)

Advance solver by one step.

#### Usage

    Lorenz_Solver$step()

------------------------------------------------------------------------

### Method `reset()`

Reset solver to its initial state.

#### Usage

    Lorenz_Solver$reset()

------------------------------------------------------------------------

### Method `collect()`

Get or set history collection behavior.

#### Usage

    Lorenz_Solver$collect(value)

#### Arguments

- `value`:

  Logical flag for collect-history behavior.

------------------------------------------------------------------------

### Method `history_size()`

Return number of stored history entries.

#### Usage

    Lorenz_Solver$history_size()

------------------------------------------------------------------------

### Method `history_step()`

Return one history entry by index.

#### Usage

    Lorenz_Solver$history_step(i)

#### Arguments

- `i`:

  Integer index into solver history.

------------------------------------------------------------------------

### Method [`history()`](https://rdrr.io/r/utils/savehistory.html)

Return history as a tibble.

#### Usage

    Lorenz_Solver$history()

------------------------------------------------------------------------

### Method `set_target()`

Set calibration targets for fitting.

#### Usage

    Lorenz_Solver$set_target(times, target, obs_indices)

#### Arguments

- `times`:

  Numeric vector of requested output times.

- `target`:

  Numeric matrix of target observations.

- `obs_indices`:

  Integer vector of observed-state indices.

------------------------------------------------------------------------

### Method `fit()`

Fit initial conditions and/or parameters.

#### Usage

    Lorenz_Solver$fit(ic = NULL, params = NULL)

#### Arguments

- `ic`:

  Optional initial condition value(s) for fitting.

- `params`:

  Optional parameter vector for fitting.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Lorenz_Solver$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
