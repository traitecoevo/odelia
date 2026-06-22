.onLoad <- function(libname, pkgname) {
  # Re-load odelia's shared library with global symbol visibility so that
  # packages that `LinkingTo: odelia` (and instantiate odelia's `Solver`) can
  # resolve the compiled XAD runtime symbols (e.g. `xad::Tape`'s destructor),
  # which live in odelia.so rather than the headers.
  #
  # R loads an imported namespace and runs its `.onLoad` before it loads the
  # importing package's own DLL, so doing this here fixes every consumer for
  # both `library()` and `devtools::load_all()` without per-package, OS-specific
  # linker hacks. See ARCHITECTURE.md for the full per-platform contract (and the
  # open Windows case, #29). Do not remove this global load without coordinating
  # with downstream consumers (e.g. plant).
  odelia_load_dll(local = FALSE)
}
