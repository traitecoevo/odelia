# What a solve cost

The number of right-hand-side evaluations (\`n_rhs\`), Jacobian
formations (\`n_jac\`), accepted steps (\`n_steps\`) and rejected
attempts (\`n_rejections\`) behind a result of \[ode_solve()\], as a
named vector. For a right-hand side that is itself expensive these are
the whole cost of the integration; \`n_rejections\` is what the
step-size controller wasted.

## Usage

``` r
ode_counts(x, ...)

# S3 method for class 'odelia_solution'
ode_counts(x, ...)
```

## Arguments

- x:

  A result of \[ode_solve()\].

- ...:

  Ignored.

## Value

A named numeric vector.

## Examples

``` r
out <- ode_solve(function(t, y, p) -y, y0 = 1, times = c(0, 1))
ode_counts(out)
#>        n_rhs        n_jac      n_steps n_rejections 
#>           73            0           12            0 
```
