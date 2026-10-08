# ODE Solver with Automatic Differentiation, in C++ Header Files

odelia provides an ODE solver in C++ header files, offering both an
explicit adaptive-step Runge-Kutta 4-5 method and an implicit,
adaptive-step RODAS4(3) Rosenbrock method for stiff systems, with an
interface to R via Rcpp. ODE systems can be templated on their scalar
type to support automatic differentiation, enabling exact reverse-mode
gradients of a recorded run and exact Jacobians for the implicit solver.
A system written as an R function can be solved by the same steppers.

## Details

The DESCRIPTION file: This package was not yet installed at build
time.  

## Package Content

Index: This package was not yet installed at build time.  

## Author

Daniel Falster \[aut, cre\] (ORCID:
\<https://orcid.org/0000-0002-9814-092X\>), Richard FitzJohn \[aut\],
Andrew O'Reilly-Nugent \[aut\]

## Maintainer

Daniel Falster \<daniel.falster@unsw.edu.au\>
