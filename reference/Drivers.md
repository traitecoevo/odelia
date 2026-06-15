# Drivers R6 Class

R6 wrapper for external forcing data (drivers)

## Public fields

- `ptr`:

  External pointer to the underlying C++ Drivers object.

## Methods

### Public methods

- [`Drivers$new()`](#method-Drivers-new)

- [`Drivers$set_constant()`](#method-Drivers-set_constant)

- [`Drivers$set_variable()`](#method-Drivers-set_variable)

- [`Drivers$set_extrapolate()`](#method-Drivers-set_extrapolate)

- [`Drivers$evaluate()`](#method-Drivers-evaluate)

- [`Drivers$evaluate_range()`](#method-Drivers-evaluate_range)

- [`Drivers$get_names()`](#method-Drivers-get_names)

- [`Drivers$clear()`](#method-Drivers-clear)

- [`Drivers$clone()`](#method-Drivers-clone)

------------------------------------------------------------------------

### Method `new()`

Initialize a \`Drivers\` instance.

#### Usage

    Drivers$new()

------------------------------------------------------------------------

### Method `set_constant()`

Set a constant driver value.

#### Usage

    Drivers$set_constant(name, value)

#### Arguments

- `name`:

  Driver variable name.

- `value`:

  Constant driver value.

------------------------------------------------------------------------

### Method `set_variable()`

Set a variable driver from paired vectors.

#### Usage

    Drivers$set_variable(name, x, y)

#### Arguments

- `name`:

  Driver variable name.

- `x`:

  Numeric vector of input coordinates (e.g., time).

- `y`:

  Numeric vector of values corresponding to \`x\`.

------------------------------------------------------------------------

### Method `set_extrapolate()`

Configure extrapolation behavior for a driver.

#### Usage

    Drivers$set_extrapolate(name, extrapolate)

#### Arguments

- `name`:

  Driver variable name.

- `extrapolate`:

  Logical flag controlling extrapolation beyond data range.

------------------------------------------------------------------------

### Method `evaluate()`

Evaluate a driver at one input value.

#### Usage

    Drivers$evaluate(name, x)

#### Arguments

- `name`:

  Driver variable name.

- `x`:

  Numeric vector of input coordinates (e.g., time).

------------------------------------------------------------------------

### Method `evaluate_range()`

Evaluate a driver over a vector of input values.

#### Usage

    Drivers$evaluate_range(name, x)

#### Arguments

- `name`:

  Driver variable name.

- `x`:

  Numeric vector of input coordinates (e.g., time).

------------------------------------------------------------------------

### Method `get_names()`

Return names of configured drivers.

#### Usage

    Drivers$get_names()

------------------------------------------------------------------------

### Method `clear()`

Remove all configured drivers.

#### Usage

    Drivers$clear()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Drivers$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
