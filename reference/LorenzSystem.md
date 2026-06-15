# Lorenz System R6 Class

R6 wrapper for Lorenz system

## Public fields

- `ptr`:

  External pointer to the underlying C++ Lorenz system object.

## Methods

### Public methods

- [`LorenzSystem$new()`](#method-LorenzSystem-new)

- [`LorenzSystem$pars()`](#method-LorenzSystem-pars)

- [`LorenzSystem$set_params()`](#method-LorenzSystem-set_params)

- [`LorenzSystem$set_state()`](#method-LorenzSystem-set_state)

- [`LorenzSystem$state()`](#method-LorenzSystem-state)

- [`LorenzSystem$rates()`](#method-LorenzSystem-rates)

- [`LorenzSystem$set_initial_state()`](#method-LorenzSystem-set_initial_state)

- [`LorenzSystem$reset()`](#method-LorenzSystem-reset)

- [`LorenzSystem$clone()`](#method-LorenzSystem-clone)

------------------------------------------------------------------------

### Method `new()`

Initialize a Lorenz system object.

#### Usage

    LorenzSystem$new(sigma, R, b)

#### Arguments

- `sigma`:

  Lorenz parameter sigma.

- `R`:

  Lorenz parameter R.

- `b`:

  Lorenz parameter b.

------------------------------------------------------------------------

### Method `pars()`

Return current system parameters.

#### Usage

    LorenzSystem$pars()

------------------------------------------------------------------------

### Method `set_params()`

Set model parameters.

#### Usage

    LorenzSystem$set_params(params)

#### Arguments

- `params`:

  Numeric vector of system parameters.

------------------------------------------------------------------------

### Method `set_state()`

Set system state and time.

#### Usage

    LorenzSystem$set_state(y, time = 0)

#### Arguments

- `y`:

  Numeric state vector.

- `time`:

  Scalar time value.

------------------------------------------------------------------------

### Method `state()`

Return current state.

#### Usage

    LorenzSystem$state()

------------------------------------------------------------------------

### Method `rates()`

Return current rates.

#### Usage

    LorenzSystem$rates()

------------------------------------------------------------------------

### Method `set_initial_state()`

Set initial state and initial time.

#### Usage

    LorenzSystem$set_initial_state(y, t0 = 0)

#### Arguments

- `y`:

  Numeric state vector.

- `t0`:

  Initial time value.

------------------------------------------------------------------------

### Method `reset()`

Reset the system to its initial condition.

#### Usage

    LorenzSystem$reset()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    LorenzSystem$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
