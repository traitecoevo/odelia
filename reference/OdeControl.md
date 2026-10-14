# OdeControl R6 Class

R6 wrapper for ODE solver control settings

## Public fields

- `ptr`:

  External pointer to the underlying C++ OdeControl object.

## Methods

### Public methods

- [`OdeControl$new()`](#method-OdeControl-initialize)

- [`OdeControl$get_controls()`](#method-OdeControl-get_controls)

- [`OdeControl$set_controls()`](#method-OdeControl-set_controls)

- [`OdeControl$set_tol_abs()`](#method-OdeControl-set_tol_abs)

- [`OdeControl$set_tol_rel()`](#method-OdeControl-set_tol_rel)

- [`OdeControl$set_a_y()`](#method-OdeControl-set_a_y)

- [`OdeControl$set_a_dydt()`](#method-OdeControl-set_a_dydt)

- [`OdeControl$set_step_size_min()`](#method-OdeControl-set_step_size_min)

- [`OdeControl$set_step_size_max()`](#method-OdeControl-set_step_size_max)

- [`OdeControl$set_controller()`](#method-OdeControl-set_controller)

- [`OdeControl$get_controller()`](#method-OdeControl-get_controller)

- [`OdeControl$set_stabilisation()`](#method-OdeControl-set_stabilisation)

- [`OdeControl$get_stabilisation()`](#method-OdeControl-get_stabilisation)

- [`OdeControl$set_step_size_initial()`](#method-OdeControl-set_step_size_initial)

- [`OdeControl$clone()`](#method-OdeControl-clone)

------------------------------------------------------------------------

### `OdeControl$new()`

Initialize an \`OdeControl\` instance.

#### Usage

    OdeControl$new()

------------------------------------------------------------------------

### `OdeControl$get_controls()`

Get current solver control settings.

#### Usage

    OdeControl$get_controls()

------------------------------------------------------------------------

### `OdeControl$set_controls()`

Set all solver control settings at once.

#### Usage

    OdeControl$set_controls(
      tol_abs,
      tol_rel,
      a_y,
      a_dydt,
      step_size_min,
      step_size_max,
      step_size_initial
    )

#### Arguments

- `tol_abs`:

  Absolute tolerance.

- `tol_rel`:

  Relative tolerance.

- `a_y`:

  State scaling coefficient.

- `a_dydt`:

  Derivative scaling coefficient.

- `step_size_min`:

  Minimum allowed step size.

- `step_size_max`:

  Maximum allowed step size.

- `step_size_initial`:

  Initial step size.

------------------------------------------------------------------------

### `OdeControl$set_tol_abs()`

Set absolute tolerance.

#### Usage

    OdeControl$set_tol_abs(tol_abs)

#### Arguments

- `tol_abs`:

  Absolute tolerance.

------------------------------------------------------------------------

### `OdeControl$set_tol_rel()`

Set relative tolerance.

#### Usage

    OdeControl$set_tol_rel(tol_rel)

#### Arguments

- `tol_rel`:

  Relative tolerance.

------------------------------------------------------------------------

### `OdeControl$set_a_y()`

Set state scaling coefficient.

#### Usage

    OdeControl$set_a_y(a_y)

#### Arguments

- `a_y`:

  State scaling coefficient.

------------------------------------------------------------------------

### `OdeControl$set_a_dydt()`

Set derivative scaling coefficient.

#### Usage

    OdeControl$set_a_dydt(a_dydt)

#### Arguments

- `a_dydt`:

  Derivative scaling coefficient.

------------------------------------------------------------------------

### `OdeControl$set_step_size_min()`

Set minimum step size.

#### Usage

    OdeControl$set_step_size_min(step_size_min)

#### Arguments

- `step_size_min`:

  Minimum allowed step size.

------------------------------------------------------------------------

### `OdeControl$set_step_size_max()`

Set maximum step size.

#### Usage

    OdeControl$set_step_size_max(step_size_max)

#### Arguments

- `step_size_max`:

  Maximum allowed step size.

------------------------------------------------------------------------

### `OdeControl$set_controller()`

Set the step-size rule, \`"gsl"\` or \`"hairer"\`.

#### Usage

    OdeControl$set_controller(controller)

#### Arguments

- `controller`:

  Which rule turns the error estimate into the next step: \`"hairer"\`
  (the default; the error in the RMS norm over components, accept at 1
  and rescale after every step, as Hairer's codes and \`deSolve\` do,
  damped by the previous step's error as \`dopri5.f\` is and
  \`deSolve\`'s \`rk\` is not; see \`set_stabilisation()\`) or \`"gsl"\`
  (the worst component's error, accept below 1.1, grow only below 0.5,
  otherwise keep the step). \`"gsl"\` rejects many more attempts;
  \`"hairer"\`'s RMS norm bounds each component more loosely when the
  state is large and the error sits in a few components.

------------------------------------------------------------------------

### `OdeControl$get_controller()`

The step-size rule in force.

#### Usage

    OdeControl$get_controller()

------------------------------------------------------------------------

### `OdeControl$set_stabilisation()`

Set the \`"hairer"\` rule's stabilisation, \`beta\` in \`dopri5.f\`:
0.04 (the default) damps step growth while the error is rising, which
keeps an explicit stepper off its stability limit; 0 is the classical
rule, which takes fewer steps on smooth problems.

#### Usage

    OdeControl$set_stabilisation(stabilisation)

#### Arguments

- `stabilisation`:

  A number in \`\[0, 0.2)\`.

------------------------------------------------------------------------

### `OdeControl$get_stabilisation()`

The \`"hairer"\` rule's stabilisation in force.

#### Usage

    OdeControl$get_stabilisation()

------------------------------------------------------------------------

### `OdeControl$set_step_size_initial()`

Set initial step size.

#### Usage

    OdeControl$set_step_size_initial(step_size_initial)

#### Arguments

- `step_size_initial`:

  Initial step size.

------------------------------------------------------------------------

### `OdeControl$clone()`

The objects of this class are cloneable with this method.

#### Usage

    OdeControl$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
