## Odelia 0.0.0.9000

Odelia is a new package, arsising out of https://github.com/traitecoevo/plant/. In that project, Rich FitzJohn built a custom ODE solver, using Runge-Kutta4-5 method, in C++. I'm spinning that code out into a package, as I want to use it elsewhere. 

* New implicit, adaptive-step **RODAS4(3)** Rosenbrock stepper for stiff systems, selectable via `method = "rodas"` when constructing a solver (#35). It reuses the existing adaptive step-size controller and obtains an exact Jacobian by forward-mode automatic differentiation; systems opt in by providing a `template<class U> rebind<U> rebind_from()` method (the same double->AD lift the gradient driver uses). The explicit RKCK 4(5) method (`method = "rkck"`) remains the default.

* `odelia` now loads its shared library with global symbol visibility in `.onLoad`, so packages that `LinkingTo: odelia` and instantiate `Solver` can resolve the compiled XAD runtime symbols at load time without per-package linker hacks (#26).

