# Signal that a state is outside the model's domain

Call this from a right-hand side, Jacobian or validity predicate handed
to \[OdeSolver\] or \[ode_solve()\], at any depth below it, to say that
the state it was given is one the model has no meaning for – a pool
below zero, a probability above one, an inner solve that did not
converge – as distinct from an error in the code. The callback is left
at once, as by \`return()\` (its \`on.exit()\` expressions run), and the
solver \*rejects the current step\* and retries it smaller, exactly as
it does for a compiled system throwing \`util::DomainError\`. Any other
error raised in a callback propagates unchanged and ends the solve,
which is what a bug should do.

## Usage

``` r
domain_error(message)
```

## Arguments

- message:

  What is wrong with the state.

## Value

Does not return.

## Details

If the smallest permitted step still leaves the domain, the solve stops
with an error naming this message, so make it say what was wrong. Called
outside a solver callback, this is an ordinary error of class
\`odelia_domain_error\`.

## Examples

``` r
rhs <- function(t, y) {
  if (y[1] < 0) domain_error("y must stay non-negative")
  -sqrt(y[1])
}
```
