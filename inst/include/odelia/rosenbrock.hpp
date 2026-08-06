// -*-c++-*-
#ifndef ODELIA_ROSENBROCK_HPP_
#define ODELIA_ROSENBROCK_HPP_

#include <vector>
#include <cmath>
#include <cstddef>
#include <XAD/XAD.hpp>
#include <odelia/ode_util.hpp>

// ROS34PW2: a four-stage, third-order, L-stable Rosenbrock-W step (Rang &
// Angermann, 2005). A W-method tolerates an approximate Jacobian, so the
// Jacobian here is a forward difference and the same one is reused across the
// step's stages -- one factorisation per step. Intended for the gentle residual
// of an operator split, where the stiff part has already been removed, so a small
// dense solve on the fast block is all that is needed.

namespace odelia {
namespace ode {

struct ROS34PW2 {
  static constexpr int s = 4;
  static constexpr double gamma = 4.3586652150845900e-01;
  // strictly-lower stage weights a_ij and lower coupling gamma_ij (diagonal = gamma)
  static constexpr double A[4][4] = {
    {0, 0, 0, 0},
    {8.7173304301691801e-01, 0, 0, 0},
    {8.4457060015369423e-01, -1.1299064236484185e-01, 0, 0},
    {0, 0, 1.0, 0}};
  static constexpr double G[4][4] = {
    {4.3586652150845900e-01, 0, 0, 0},
    {-8.7173304301691801e-01, 4.3586652150845900e-01, 0, 0},
    {-9.0338057013044082e-01, 5.4180672388095326e-02, 4.3586652150845900e-01, 0},
    {2.4212380706095346e-01, -1.2232505839045147e+00, 5.4526025533510214e-01, 4.3586652150845900e-01}};
  static constexpr double b[4] =
    {2.4212380706095346e-01, -1.2232505839045147e+00, 1.5452602553351020e+00, 4.3586652150845900e-01};
  static constexpr double b2[4] =
    {3.7810903145819369e-01, -9.6042292212423178e-02, 5.0000000000000000e-01, 2.1793326075422950e-01};
};

// Solve M z = rhs in place by Gaussian elimination with partial pivoting; M is a
// double n x n row-major matrix (overwritten), rhs any scalar T (double or an
// active AD type). Pivots on the passive matrix, so the same factorisation drives
// an active back-substitution -- the discrete adjoint of a solve with a constant
// matrix. For the small fast block a dense solve beats any sparsity bookkeeping.
template <class T>
void dense_solve(std::vector<double>& M, std::vector<T>& rhs, int n) {
  for (int col = 0; col < n; ++col) {
    int piv = col;
    for (int r = col + 1; r < n; ++r)
      if (std::fabs(M[r * n + col]) > std::fabs(M[piv * n + col])) piv = r;
    if (piv != col) {
      for (int c = 0; c < n; ++c) std::swap(M[piv * n + c], M[col * n + c]);
      std::swap(rhs[piv], rhs[col]);
    }
    const double d = M[col * n + col];
    for (int r = col + 1; r < n; ++r) {
      const double f = M[r * n + col] / d;
      for (int c = col; c < n; ++c) M[r * n + c] -= f * M[col * n + c];
      rhs[r] -= f * rhs[col];
    }
  }
  for (int r = n - 1; r >= 0; --r) {
    T acc = rhs[r];
    for (int c = r + 1; c < n; ++c) acc -= M[r * n + c] * rhs[c];
    rhs[r] = acc / M[r * n + r];
  }
}

// One ROS34PW2 step of size h on the autonomous system y' = res(y). Updates y and
// fills yerr with the embedded (order-2) error estimate. `res` is a functor with a
// templated call operator res(const std::vector<U>&, std::vector<U>&) usable at
// both double and the active scalar S. The Jacobian is a forward difference on the
// values (a W-method tolerates an approximate one) and reused across the stages;
// held off the tape, it makes each stage a solve with a constant matrix, so the
// step differentiates by an active back-substitution with recorded factors.
template <class S, class Residual>
void ros34pw2_step(const Residual& res, std::vector<S>& y, double h, std::vector<S>& yerr) {
  using R = ROS34PW2;
  const int n = (int)y.size();

  std::vector<double> yv(n), f0(n), fp(n), yd(n);
  for (int r = 0; r < n; ++r) yv[r] = xad::value(y[r]);
  res(yv, f0);
  std::vector<double> J((size_t)n * n);
  for (int col = 0; col < n; ++col) {
    const double d = 1e-7 * std::max(1.0, std::fabs(yv[col]));
    yd = yv; yd[col] += d;
    res(yd, fp);
    for (int r = 0; r < n; ++r) J[r * n + col] = (fp[r] - f0[r]) / d;
  }
  std::vector<double> Abase((size_t)n * n);
  for (int r = 0; r < n; ++r)
    for (int c = 0; c < n; ++c)
      Abase[r * n + c] = (r == c ? 1.0 : 0.0) - h * R::gamma * J[r * n + c];

  std::vector<std::vector<S>> k(R::s, std::vector<S>(n));
  std::vector<S> ytmp(n), stagef(n), Jsum(n), rhs(n);
  for (int i = 0; i < R::s; ++i) {
    for (int r = 0; r < n; ++r) {
      ytmp[r] = y[r];
      for (int j = 0; j < i; ++j) ytmp[r] += R::A[i][j] * k[j][r];
    }
    res(ytmp, stagef);
    for (int r = 0; r < n; ++r) {
      S acc(0.0);
      for (int j = 0; j < i; ++j) {
        const double g = R::G[i][j];
        if (g != 0.0) for (int c = 0; c < n; ++c) acc += J[r * n + c] * (g * k[j][c]);
      }
      Jsum[r] = acc;
    }
    for (int r = 0; r < n; ++r) rhs[r] = h * stagef[r] + h * Jsum[r];
    std::vector<double> M = Abase;
    dense_solve(M, rhs, n);
    k[i] = rhs;
  }

  for (int r = 0; r < n; ++r) {
    S yn(0.0), err(0.0);
    for (int i = 0; i < R::s; ++i) { yn += R::b[i] * k[i][r]; err += (R::b[i] - R::b2[i]) * k[i][r]; }
    y[r] += yn;
    yerr[r] = err;
  }
}

}
}

#endif
