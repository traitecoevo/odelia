## Odelia 0.0.0.9000

Odelia is a new package, arsising out of https://github.com/traitecoevo/plant/. In that project, Rich FitzJohn built a custom ODE solver, using Runge-Kutta4-5 method, in C++. I'm spinning that code out into a package, as I want to use it elsewhere. 

* `odelia` now loads its shared library with global symbol visibility in `.onLoad`, so packages that `LinkingTo: odelia` and instantiate `Solver` can resolve the compiled XAD runtime symbols at load time without per-package linker hacks (#26).

