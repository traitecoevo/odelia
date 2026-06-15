# Load odelia shared library with configurable symbol visibility

Locates the installed \`odelia\` shared library and calls
\[base::dyn.load()\] on it. This is useful in \`Rcpp::sourceCpp()\`
workflows where symbols from the package shared object must be available
to a temporary shared object.

## Usage

``` r
odelia_load_dll(local = FALSE, now = TRUE)
```

## Arguments

- local:

  Passed to \[base::dyn.load()\]. Default \`FALSE\` to expose symbols
  globally for downstream dynamic linking.

- now:

  Passed to \[base::dyn.load()\]. Default \`TRUE\`.

## Value

Invisibly returns the resolved path to the \`odelia\` shared library.
