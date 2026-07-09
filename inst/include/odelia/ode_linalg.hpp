// -*-c++-*-
#ifndef ODELIA_ODE_LINALG_HPP_
#define ODELIA_ODE_LINALG_HPP_

// Minimal dense linear algebra for the implicit (Rosenbrock) stepper.
//
// A hand-rolled LU with partial pivoting, templated on the scalar type so it
// works unchanged for `double` and for XAD active types (AReal/FReal). The
// arithmetic is plain operator use, so under an active scalar every operation is
// recorded on the tape and gradients flow through the solve. Pivot *selection*
// compares magnitudes of the underlying numeric values via `xad::value(...)`, so
// the choice of pivot never itself becomes a taped/differentiated quantity.
//
// Matrices are stored row-major in a flat std::vector: A(i, j) == a[i * n + j].
// One factorization is reused across all six Rosenbrock stage solves.

#include <vector>
#include <cstddef>
#include <cmath>
#include <XAD/XAD.hpp>
#include <odelia/ode_util.hpp>

namespace odelia {
namespace ode {
namespace linalg {

// In-place LU decomposition with partial pivoting (Doolittle, unit lower
// diagonal). On return `a` holds L (below the diagonal, implicit unit diagonal)
// and U (on and above the diagonal), and `piv` holds the row permutation such
// that row i of the factorization corresponds to original row piv[i].
template <typename T>
void lu_decompose(std::vector<T>& a, size_t n, std::vector<size_t>& piv) {
  piv.resize(n);
  for (size_t i = 0; i < n; ++i) {
    piv[i] = i;
  }

  for (size_t k = 0; k < n; ++k) {
    // Find pivot row: largest |a(i,k)| for i >= k, comparing numeric values.
    size_t pivot_row = k;
    double best = std::abs(xad::value(a[k * n + k]));
    for (size_t i = k + 1; i < n; ++i) {
      const double v = std::abs(xad::value(a[i * n + k]));
      if (v > best) {
        best = v;
        pivot_row = i;
      }
    }

    if (best == 0.0) {
      util::stop("Singular matrix in LU decomposition (implicit stepper)");
    }

    // Swap rows k and pivot_row (both in the matrix and the permutation).
    if (pivot_row != k) {
      for (size_t j = 0; j < n; ++j) {
        std::swap(a[k * n + j], a[pivot_row * n + j]);
      }
      std::swap(piv[k], piv[pivot_row]);
    }

    // Eliminate below the pivot.
    const T inv_pivot = T(1.0) / a[k * n + k];
    for (size_t i = k + 1; i < n; ++i) {
      const T factor = a[i * n + k] * inv_pivot;
      a[i * n + k] = factor; // store multiplier (L)
      for (size_t j = k + 1; j < n; ++j) {
        a[i * n + j] -= factor * a[k * n + j];
      }
    }
  }
}

// Solve A x = b given the LU factorization from lu_decompose. `b` is the
// right-hand side (unpermuted); the solution is written into `x`.
template <typename T>
void lu_solve(const std::vector<T>& lu, size_t n,
              const std::vector<size_t>& piv,
              const std::vector<T>& b, std::vector<T>& x) {
  x.resize(n);

  // Forward substitution (Ly = Pb), applying the row permutation to b.
  for (size_t i = 0; i < n; ++i) {
    T sum = b[piv[i]];
    for (size_t j = 0; j < i; ++j) {
      sum -= lu[i * n + j] * x[j];
    }
    x[i] = sum; // unit lower diagonal, no division
  }

  // Back substitution (Ux = y).
  for (size_t i = n; i-- > 0;) {
    T sum = x[i];
    for (size_t j = i + 1; j < n; ++j) {
      sum -= lu[i * n + j] * x[j];
    }
    x[i] = sum / lu[i * n + i];
  }
}

} // namespace linalg
} // namespace ode
} // namespace odelia

#endif
