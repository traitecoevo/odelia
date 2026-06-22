# odelia architecture: the XAD `Tape` runtime and how downstream packages link it

> **Read this before changing how XAD is compiled, how `odelia.so`/`odelia.dll`
> is loaded, or `R/zzz.R`.** These pieces form a contract that packages which
> `LinkingTo: odelia` depend on. Breaking it produces *link-time* or *load-time*
> failures in downstream packages (e.g. plant) that do not show up in odelia's
> own checks.

## Why there is compiled code at all

odelia is *almost* header-only — the ODE `Solver`, interpolator, and ODE
interface are templates that compile into whichever package includes them. The
exception is the **XAD automatic-differentiation `Tape` runtime**:

- `inst/include/XAD/Tape.hpp` **declares** the `Tape<T,N>` members and, near the
  bottom, does `extern template class Tape<double>;` (and friends). `extern
  template` tells every translation unit *"do not instantiate `Tape` here — it
  is compiled elsewhere."*
- `src/Tape.cpp` holds the **only** copies of the `Tape<T,N>` member-function
  *definitions* (`~Tape()`, `computeAdjoints()`, `resetTo()`, …) **and** the
  explicit instantiations (`template class Tape<double>;`, via the
  `MAKE_TAPE_TLS` macros). These compile into `odelia.so` / `odelia.dll` and
  nowhere else.

This split is a deliberate compile-time optimisation (the heavy `Tape` code is
compiled once, in odelia, not in every consumer). Its consequence is that **any
consumer that instantiates `Solver<…>` emits undefined references to
`xad::Tape<double,1>::~Tape()` et al. that must be resolved against odelia's
compiled library.**

## The single-`active_tape_` invariant

`Tape.hpp` declares `static XAD_THREAD_LOCAL Tape* active_tape_;`, defined once
in `Tape.cpp`. This is **one active-tape pointer per thread, per process** — AD
recording (`getActive()`, `activate()`, `deactivate()`) all key off it. The
whole design assumes a *single* definition of this symbol in the process.

Keeping `Tape` compiled in exactly one place (odelia's library) is what
guarantees that single definition. Any change that causes `Tape` to be
instantiated in *each* consumer (e.g. making it header-only, or shipping the
definitions as an includable impl header) gives every consumer DLL its own
`active_tape_`, which is only safe if a given AD computation never crosses the
odelia↔consumer DLL boundary — a subtle, easy-to-violate condition. Do not do
this without understanding the implication.

## How downstream packages resolve the symbols, per platform

A package that `LinkingTo: odelia` and instantiates `Solver` must get odelia's
compiled `Tape` symbols resolved. Two things make that work today:

1. **odelia loads its DLL with global symbol visibility** in `R/zzz.R`
   (`.onLoad` → `odelia_load_dll(local = FALSE)`; see odelia #26/#28). R runs an
   imported namespace's `.onLoad` *before* it loads the importing package's own
   DLL, so odelia's symbols are present in the process first.
2. **The consumer must load odelia's namespace before its own DLL.** `LinkingTo`
   is C++-only and does **not** load odelia's namespace; the consumer therefore
   needs `Imports: odelia` *and* an actual NAMESPACE import so R loads odelia
   (running the `.onLoad` above) ahead of the consumer's `useDynLib`. plant does
   this with `@importFrom odelia odelia_load_dll` (mirroring the `@importFrom
   Rcpp evalCpp` idiom). `devtools::load_all()` happens to load `Imports` first,
   but a fresh `library()` only does so via a real NAMESPACE import.

| Platform | Linker behaviour for the undefined `Tape` symbol | Status |
|---|---|---|
| **macOS** | `-undefined dynamic_lookup` (flat namespace): allowed in `.so`, bound lazily at first call; resolved against globally-loaded odelia | ✅ works (steps 1–2) |
| **Linux** | GNU `ld -shared` allows undefined symbols in `.so`, resolved at `dlopen` against globally-loaded odelia | ✅ works (steps 1–2) |
| **Windows** | mingw `ld` **requires every symbol resolved at link time**; a DLL cannot carry undefined imports, so the runtime trick in step 1 cannot help | ❌ **open — see [traitecoevo/odelia#29](https://github.com/traitecoevo/odelia/issues/29)** |

On Windows the consumer's `plant.dll` fails to link with `undefined reference to
xad::Tape<double, 1ull>::~Tape()`. The fix (link-time access to the single
compiled `Tape`, e.g. an import library or a documented link line) is tracked in
**#29**; update this table and the downstream contract below once it lands.

## Contract for `LinkingTo: odelia` consumers

- Add `odelia` to `LinkingTo:` **and** `Imports:` in `DESCRIPTION`.
- Add a real NAMESPACE import from odelia (e.g. `@importFrom odelia
  odelia_load_dll`) so odelia's namespace — and its global-load `.onLoad` —
  runs before your package's `useDynLib`.
- macOS / Linux: nothing else is required.
- Windows: pending #29.

## Do NOT change without coordinating downstream

- The `extern template` / `Tape.cpp` split, or making `Tape` header-only — see
  the single-`active_tape_` invariant above.
- The set of explicit `Tape` instantiations in `src/Tape.cpp` (consumers rely on
  the exact instantiated types).
- The global-visibility load in `R/zzz.R` (`odelia_load_dll(local = FALSE)`).

Any of these will break downstream linking/loading silently (odelia's own checks
will still pass). Coordinate with consumers (plant) and update this document.
