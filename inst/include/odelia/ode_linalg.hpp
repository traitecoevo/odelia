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
#include <algorithm>
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

// ---------------------------------------------------------------------------
// Eigenvalues of a general real matrix.
//
// Used by the steady-state solver to test whether a fixed point is attracting
// (all eigenvalues of the state Jacobian df/dy have negative real part). Only
// the *numeric* eigenvalues are needed -- the matrix is already at a converged
// equilibrium and the stability verdict is a diagnostic, not a differentiated
// quantity -- so the routine strips any active AD layer via xad::value and works
// entirely in double. This deliberately does not tape: stability is not a
// gradient target.
//
// Method: Householder reduction to upper Hessenberg form, then the Francis
// double-shift QR algorithm (eigenvalues only, no vectors). Standard dense
// eigenvalue machinery (Golub & Van Loan, Matrix Computations, ch. 7); no
// external linear-algebra library is required, matching the hand-rolled LU.
namespace detail {

// Reduce a (row-major, n*n) real matrix to upper Hessenberg form in place by
// Householder similarity transforms (eigenvalues preserved).
inline void to_hessenberg(std::vector<double>& a, size_t n) {
  auto A = [&](size_t i, size_t j) -> double& { return a[i * n + j]; };
  if (n < 3) {
    return;
  }
  std::vector<double> v(n, 0.0);
  for (size_t k = 0; k + 2 < n; ++k) {
    double scale = 0.0;
    for (size_t i = k + 1; i < n; ++i) {
      scale += std::fabs(A(i, k));
    }
    if (scale == 0.0) {
      continue; // column already reduced
    }
    double h = 0.0;
    for (size_t i = k + 1; i < n; ++i) {
      v[i] = A(i, k) / scale;
      h += v[i] * v[i];
    }
    const double f = v[k + 1];
    const double g = (f >= 0.0) ? -std::sqrt(h) : std::sqrt(h);
    h -= f * g;
    v[k + 1] = f - g;
    // Apply P = I - v v^T / h on the left (P A) then on the right (A P).
    for (size_t j = 0; j < n; ++j) {
      double s = 0.0;
      for (size_t i = k + 1; i < n; ++i) {
        s += v[i] * A(i, j);
      }
      s /= h;
      for (size_t i = k + 1; i < n; ++i) {
        A(i, j) -= s * v[i];
      }
    }
    for (size_t i = 0; i < n; ++i) {
      double s = 0.0;
      for (size_t j = k + 1; j < n; ++j) {
        s += A(i, j) * v[j];
      }
      s /= h;
      for (size_t j = k + 1; j < n; ++j) {
        A(i, j) -= s * v[j];
      }
    }
    A(k + 1, k) = g * scale;
    for (size_t i = k + 2; i < n; ++i) {
      A(i, k) = 0.0;
    }
  }
}

// Francis double-shift QR on an upper Hessenberg matrix (destroyed in place).
// Writes the real and imaginary parts of the eigenvalues to wr/wi. Returns
// false if an eigenvalue fails to converge within the iteration budget.
inline bool hqr(std::vector<double>& a, size_t n, std::vector<double>& wr,
                std::vector<double>& wi) {
  auto A = [&](int i, int j) -> double& {
    return a[static_cast<size_t>(i) * n + static_cast<size_t>(j)];
  };
  wr.assign(n, 0.0);
  wi.assign(n, 0.0);
  const double eps = 2.220446049250313e-16;

  double anorm = 0.0;
  for (int i = 0; i < static_cast<int>(n); ++i) {
    for (int j = std::max(i - 1, 0); j < static_cast<int>(n); ++j) {
      anorm += std::fabs(A(i, j));
    }
  }

  int nn = static_cast<int>(n) - 1;
  double t = 0.0;
  while (nn >= 0) {
    int its = 0;
    int l;
    do {
      // Locate a negligible subdiagonal element to split off a subproblem.
      for (l = nn; l >= 1; --l) {
        double s = std::fabs(A(l - 1, l - 1)) + std::fabs(A(l, l));
        if (s == 0.0) {
          s = anorm;
        }
        if (std::fabs(A(l, l - 1)) <= eps * s) {
          A(l, l - 1) = 0.0;
          break;
        }
      }
      double x = A(nn, nn);
      if (l == nn) { // one real eigenvalue
        wr[nn] = x + t;
        wi[nn] = 0.0;
        --nn;
      } else {
        double y = A(nn - 1, nn - 1);
        double w = A(nn, nn - 1) * A(nn - 1, nn);
        if (l == nn - 1) { // a trailing 2x2 block: two eigenvalues
          double p = 0.5 * (y - x);
          double q = p * p + w;
          double z = std::sqrt(std::fabs(q));
          x += t;
          if (q >= 0.0) { // real pair
            z = p + (p >= 0.0 ? std::fabs(z) : -std::fabs(z));
            wr[nn - 1] = wr[nn] = x + z;
            if (z != 0.0) {
              wr[nn] = x - w / z;
            }
            wi[nn - 1] = wi[nn] = 0.0;
          } else { // complex conjugate pair
            wr[nn - 1] = wr[nn] = x + p;
            wi[nn - 1] = -z;
            wi[nn] = z;
          }
          nn -= 2;
        } else { // no eigenvalue split off: do a double-shift QR sweep
          if (its >= 60) {
            return false;
          }
          if (its == 10 || its == 20 || its == 30) { // exceptional shift
            t += x;
            for (int i = 0; i <= nn; ++i) {
              A(i, i) -= x;
            }
            double s = std::fabs(A(nn, nn - 1)) + std::fabs(A(nn - 1, nn - 2));
            y = x = 0.75 * s;
            w = -0.4375 * s * s;
          }
          ++its;
          double p = 0.0, q = 0.0, r = 0.0;
          int m;
          for (m = nn - 2; m >= l; --m) {
            double z = A(m, m);
            r = x - z;
            double s = y - z;
            p = (r * s - w) / A(m + 1, m) + A(m, m + 1);
            q = A(m + 1, m + 1) - z - r - s;
            r = A(m + 2, m + 1);
            double sc = std::fabs(p) + std::fabs(q) + std::fabs(r);
            p /= sc;
            q /= sc;
            r /= sc;
            if (m == l) {
              break;
            }
            double u = std::fabs(A(m, m - 1)) * (std::fabs(q) + std::fabs(r));
            double vv = std::fabs(p) * (std::fabs(A(m - 1, m - 1)) +
                                        std::fabs(z) + std::fabs(A(m + 1, m + 1)));
            if (u <= eps * vv) {
              break;
            }
          }
          for (int i = m + 2; i <= nn; ++i) {
            A(i, i - 2) = 0.0;
            if (i != m + 2) {
              A(i, i - 3) = 0.0;
            }
          }
          for (int k = m; k <= nn - 1; ++k) {
            if (k != m) {
              p = A(k, k - 1);
              q = A(k + 1, k - 1);
              r = 0.0;
              if (k + 1 != nn) {
                r = A(k + 2, k - 1);
              }
              x = std::fabs(p) + std::fabs(q) + std::fabs(r);
              if (x != 0.0) {
                p /= x;
                q /= x;
                r /= x;
              }
            }
            double s = (p >= 0.0 ? 1.0 : -1.0) * std::sqrt(p * p + q * q + r * r);
            if (s == 0.0) {
              continue;
            }
            if (k == m) {
              if (l != m) {
                A(k, k - 1) = -A(k, k - 1);
              }
            } else {
              A(k, k - 1) = -s * x;
            }
            p += s;
            const double px = p / s, qx = q / s, rx = r / s;
            const double qp = q / p, rp = r / p;
            for (int j = k; j <= nn; ++j) { // row transform
              double pp = A(k, j) + qp * A(k + 1, j);
              if (k + 1 != nn) {
                pp += rp * A(k + 2, j);
                A(k + 2, j) -= pp * rx;
              }
              A(k + 1, j) -= pp * qx;
              A(k, j) -= pp * px;
            }
            const int mmax = (nn < k + 3) ? nn : k + 3;
            for (int i = l; i <= mmax; ++i) { // column transform
              double pp = px * A(i, k) + qx * A(i, k + 1);
              if (k + 1 != nn) {
                pp += rx * A(i, k + 2);
                A(i, k + 2) -= pp * rp;
              }
              A(i, k + 1) -= pp * qp;
              A(i, k) -= pp;
            }
          }
        }
      }
    } while (l < nn - 1);
  }
  return true;
}

} // namespace detail

// Eigenvalues of a general real n*n matrix `A` (row-major, any scalar type;
// active layers are stripped via xad::value). Real and imaginary parts are
// written to `wr` and `wi`. Throws if the QR iteration fails to converge.
template <typename T>
void eigenvalues(const std::vector<T>& A, size_t n, std::vector<double>& wr,
                 std::vector<double>& wi) {
  std::vector<double> h(n * n);
  for (size_t i = 0; i < n * n; ++i) {
    h[i] = xad::value(A[i]);
  }
  detail::to_hessenberg(h, n);
  if (!detail::hqr(h, n, wr, wi)) {
    util::stop("Eigenvalue QR iteration failed to converge");
  }
}

// The spectral abscissa: the largest real part over all eigenvalues of `A`. A
// fixed point of an autonomous system is asymptotically stable (attracting) iff
// this is < 0 for the state Jacobian df/dy.
template <typename T>
double spectral_abscissa(const std::vector<T>& A, size_t n) {
  std::vector<double> wr, wi;
  eigenvalues(A, n, wr, wi);
  double amax = wr.empty() ? 0.0 : wr[0];
  for (double re : wr) {
    if (re > amax) {
      amax = re;
    }
  }
  return amax;
}

} // namespace linalg
} // namespace ode
} // namespace odelia

#endif
