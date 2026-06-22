# Vendored XAD and local patches

odelia vendors the [XAD](https://github.com/auto-differentiation/xad) automatic
differentiation library:

| Upstream                | In this package          |
| ----------------------- | ------------------------ |
| `src/XAD/*.hpp`         | `inst/include/XAD/*.hpp` |
| `src/Tape.cpp`          | `src/Tape.cpp`           |

Some headers (`Config.hpp`, `Version.hpp`, `Instantiations.hpp`) are **generated
by XAD's CMake build**, so a faithful re-vendor is not a pure file copy — you
must configure XAD with CMake once and copy the generated headers alongside the
hand-written ones. The currently vendored snapshot is XAD `1.9.0-dev`.

## Local patches

We carry local fixes on top of the vendored snapshot, stored as `git apply`-able
patches in [`patches/`](patches/):

- **`0001-xad-libcxx21-promote-guard.patch`** — libc++ 21 (`_LIBCPP_VERSION >=
  210000`) removed the specializable `std::__promote` class template that XAD
  specializes in `Literals.hpp`, replacing it with a SFINAE-friendly
  `__promote_t` alias. XAD's specializations are unnecessary and ill-formed
  there, so the patch restricts them to `_LIBCPP_VERSION < 210000`. Not yet
  fixed upstream — drop this patch once upstream handles the rename.

## Re-vendoring workflow

1. Copy the new upstream `src/XAD/*.hpp` + CMake-generated headers into
   `inst/include/XAD/`, and `src/Tape.cpp` into `src/`.
2. Re-apply the local patches:

   ```sh
   tools/apply-xad-patches.sh
   ```

3. Rebuild and test: `make test-installed`.

Before relying on a patch, you can confirm it still applies to the freshly
vendored sources with `tools/apply-xad-patches.sh --check`. If a patch no longer
applies, upstream changed the relevant code — check whether the underlying issue
is fixed upstream (in which case delete the patch) or re-create it.
