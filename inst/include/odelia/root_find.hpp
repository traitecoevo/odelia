// -*-c++-*-
#ifndef ODELIA_ROOT_FIND_HPP_
#define ODELIA_ROOT_FIND_HPP_

// Bracketed one-dimensional root-finding on a target double -> double: TOMS748
// (Alefeld, Potra and Shi 1995, ACM TOMS 21(3), 327-344) and bisection. Each
// keeps a bracket [a, b] across which the target changes sign. TOMS748
// converges superlinearly on a smooth target; bisection gains one bit per
// evaluation whatever the target.
//
// Both are Boost.Math's (`toms748_solve` in boost/math/tools/toms748_solve.hpp
// and `bisect` in boost/math/tools/roots.hpp, Boost 1.90), carried here in
// double so that the header core needs no Boost. The arithmetic is Boost's, in
// Boost's order and split across the same functions, so a search evaluates the
// target at the same points and ends on the same bracket as Boost's does;
// tests/testthat/test-root-find.R holds the two to that. One condition: the
// target is called by reference, where Boost calls a fresh copy at every step,
// so a target keeping state of its own by value searches differently. Boost's
// code is under the Boost Software License, reproduced at the end of this
// file.
//
// What is odelia's own is the refusal. A search that ends with its bracket
// wider than the tolerance throws util::RootFindError. Boost returns that
// bracket as though it had converged, and the iteration count it reports
// cannot say whether it did: the count is the iterations used, which never
// exceed the cap.
//
// A search allocates nothing and takes the target by reference, so it can run
// on any callable, the state a step's dense output gives at a time included.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <odelia/ode_util.hpp>

namespace odelia {
namespace root {

// When a bracket [a, b] is narrow enough,
//     |a - b| < atol + rtol * min(|a|, |b|),
// and how many evaluations of the target a search may spend reaching that. The
// cap is a backstop: a target a method suits meets any usable tolerance far
// inside it.
struct tolerance {
  double atol;
  double rtol;
  std::size_t max_iterations = 1000;

  bool met(double a, double b) const {
    return std::abs(a - b) < atol + rtol * std::min(std::abs(a), std::abs(b));
  }
};

namespace internal {

// Out of line and cold: building a message inlined into a search that runs
// inside a step costs the step, in code size and in what the compiler will
// inline around it.
[[noreturn, gnu::cold, gnu::noinline]] inline void
refuse_bracket(const char* method, const char* why, double a, double b,
               double fa, double fb) {
  throw util::BracketError(
      std::string(method) + ": " + why + ": f(" + util::format_double(a) +
      ") = " + util::format_double(fa) + ", f(" + util::format_double(b) +
      ") = " + util::format_double(fb));
}

[[noreturn, gnu::cold, gnu::noinline]] inline void
refuse_unconverged(const char* method, double a, double b,
                   const tolerance& tol, std::uintmax_t used) {
  throw util::RootFindError(
      std::string(method) + ": the bracket [" + util::format_double(a) + ", " +
      util::format_double(b) + "] is still wider than atol " +
      util::format_double(tol.atol) + " + rtol " +
      util::format_double(tol.rtol) + " * min(|a|, |b|) after " +
      util::to_string(used) + " of " + util::to_string(tol.max_iterations) +
      " iterations");
}

// Ends on the root (first == second, where the target was exactly zero), or a
// bracket narrower than the tolerance: anything else is refused.
inline void check_converged(const char* method,
                            const std::pair<double, double>& r,
                            const tolerance& tol, std::uintmax_t used) {
  if (!(r.first == r.second || tol.met(r.first, r.second))) [[unlikely]] {
    refuse_unconverged(method, r.first, r.second, tol, used);
  }
}

// boost::math::sign: 0 at zero, otherwise -1 or 1 by the sign bit, which a NaN
// has too.
inline int sign(double z) {
  return (z == 0) ? 0 : std::signbit(z) ? -1 : 1;
}

// Evaluates f at c inside [a, b] and keeps the half that still brackets the
// root, or sets a = c where f(c) is exactly zero. d and fd receive the end that
// was dropped, the third-best estimate of the root. A c within a few ulps of
// either end is moved inside first. Not on a bracket of a few ulps below zero,
// as in Boost: the first test needs a > 0, and the second can then put c past
// b, which the bracket keeps.
template <class F>
void bracket(F& f, double& a, double& b, double c, double& fa, double& fb,
             double& d, double& fd) {
  double tol = std::numeric_limits<double>::epsilon() * 2;
  if ((b - a) < 2 * tol * a) {
    c = a + (b - a) / 2;
  } else if (c <= a + std::fabs(a) * tol) {
    c = a + std::fabs(a) * tol;
  } else if (c >= b - std::fabs(b) * tol) {
    c = b - std::fabs(b) * tol;
  }
  double fc = f(c);
  if (fc == 0) {
    a = c;
    fa = 0;
    d = 0;
    fd = 0;
    return;
  }
  if (sign(fa) * sign(fc) < 0) {
    d = b;
    fd = fb;
    b = c;
    fb = fc;
  } else {
    d = a;
    fd = fa;
    a = c;
    fa = fc;
  }
}

// num / denom, or r where the quotient would overflow.
inline double safe_div(double num, double denom, double r) {
  if (std::fabs(denom) < 1) {
    if (std::fabs(denom * std::numeric_limits<double>::max()) <=
        std::fabs(num)) {
      return r;
    }
  }
  return num / denom;
}

// The secant through (a, fa) and (b, fb), or the midpoint where that lands
// within a few ulps of either end: this is only reached after an
// interpolation has failed, so a root that close to an end is unlikely.
inline double secant_interpolate(const double& a, const double& b,
                                 const double& fa, const double& fb) {
  double tol = std::numeric_limits<double>::epsilon() * 5;
  double c = a - (fa / (fb - fa)) * (b - a);
  if ((c <= a + std::fabs(a) * tol) || (c >= b - std::fabs(b) * tol)) {
    return (a + b) / 2;
  }
  return c;
}

// The root of the quadratic through a, b and d (d outside [a, b]), by `count`
// Newton steps from the end where the quadratic's curvature and f agree in
// sign; a secant step where that fails or leaves (a, b).
inline double quadratic_interpolate(const double& a, const double& b,
                                    double const& d, const double& fa,
                                    const double& fb, double const& fd,
                                    unsigned count) {
  double B = safe_div(double(fb - fa), double(b - a),
                      std::numeric_limits<double>::max());
  double A = safe_div(double(fd - fb), double(d - b),
                      std::numeric_limits<double>::max());
  A = safe_div(double(A - B), double(d - a), double(0));
  if (A == 0) {
    return secant_interpolate(a, b, fa, fb);
  }
  double c;
  if (sign(A) * sign(fa) > 0) {
    c = a;
  } else {
    c = b;
  }
  for (unsigned i = 1; i <= count; ++i) {
    c -= safe_div(double(fa + (B + A * (c - b)) * (c - a)),
                  double(B + A * (2 * c - a - b)), double(1 + c - a));
  }
  if ((c <= a) || (c >= b)) {
    c = secant_interpolate(a, b, fa, fb);
  }
  return c;
}

// Inverse cubic interpolation through a, b, d and e (d and e outside [a, b]);
// quadratic interpolation where that leaves (a, b).
inline double cubic_interpolate(const double& a, const double& b,
                                const double& d, const double& e,
                                const double& fa, const double& fb,
                                const double& fd, const double& fe) {
  double q11 = (d - e) * fd / (fe - fd);
  double q21 = (b - d) * fb / (fd - fb);
  double q31 = (a - b) * fa / (fb - fa);
  double d21 = (b - d) * fd / (fd - fb);
  double d31 = (a - b) * fb / (fb - fa);
  double q22 = (d21 - q11) * fb / (fe - fb);
  double q32 = (d31 - q21) * fa / (fd - fa);
  double d32 = (d31 - q21) * fd / (fd - fa);
  double q33 = (d32 - q22) * fa / (fe - fa);
  double c = q31 + q32 + q33 + a;
  if ((c <= a) || (c >= b)) {
    c = quadratic_interpolate(a, b, d, fa, fb, fd, 3);
  }
  return c;
}

// Boost's toms748_solve, with the caller's end values. Returns the bracket it
// ended on, and leaves in `max_iter` the evaluations it spent. Refuses ends that
// do not bracket a root, but not a search that runs out: the public forms do.
// Ends meeting the tolerance are returned before their signs are read, as in
// Boost.
//
// One departure: with no evaluations to spend, Boost returns [ax, bx]
// unchecked. Here the ends are checked first, so a cap of 0 refuses a bad
// bracket and returns an exact zero at an end, as every other cap does.
template <class F>
std::pair<double, double> toms748_solve(F& f, const double& ax,
                                        const double& bx, const double& fax,
                                        const double& fbx, const tolerance& tol,
                                        std::uintmax_t& max_iter) {
  std::uintmax_t count = max_iter;
  double a, b, fa, fb, c, u, fu, a0, b0, d, fd, e, fe;
  static const double mu = 0.5f;

  a = ax;
  b = bx;
  if (a >= b) [[unlikely]] {
    refuse_bracket("toms748", "the ends are out of order", ax, bx, fax, fbx);
  }
  fa = fax;
  fb = fbx;

  if (tol.met(a, b) || (fa == 0) || (fb == 0)) {
    max_iter = 0;
    if (fa == 0) {
      b = a;
    } else if (fb == 0) {
      a = b;
    }
    return std::make_pair(a, b);
  }

  if (sign(fa) * sign(fb) > 0) [[unlikely]] {
    refuse_bracket("toms748", "the ends do not bracket a root", ax, bx, fax,
                   fbx);
  }
  if (count == 0) {
    return std::make_pair(a, b);
  }
  // Placeholders for the third and fourth estimates the first steps lack.
  fe = e = fd = 1e5F;

  if (fa != 0) {
    // A secant step first, then a quadratic one.
    c = secant_interpolate(a, b, fa, fb);
    bracket(f, a, b, c, fa, fb, d, fd);
    --count;

    if (count && (fa != 0) && !tol.met(a, b)) {
      c = quadratic_interpolate(a, b, d, fa, fb, fd, 2);
      e = d;
      fe = fd;
      bracket(f, a, b, c, fa, fb, d, fd);
      --count;
    }
  }

  while (count && (fa != 0) && !tol.met(a, b)) {
    a0 = a;
    b0 = b;
    // Cubic interpolation needs fa, fb, fd and fe distinct; where any two are
    // not, a quadratic step instead.
    double min_diff = std::numeric_limits<double>::min() * 32;
    bool prof = (std::fabs(fa - fb) < min_diff) ||
                (std::fabs(fa - fd) < min_diff) ||
                (std::fabs(fa - fe) < min_diff) ||
                (std::fabs(fb - fd) < min_diff) ||
                (std::fabs(fb - fe) < min_diff) ||
                (std::fabs(fd - fe) < min_diff);
    if (prof) {
      c = quadratic_interpolate(a, b, d, fa, fb, fd, 2);
    } else {
      c = cubic_interpolate(a, b, d, e, fa, fb, fd, fe);
    }
    e = d;
    fe = fd;
    bracket(f, a, b, c, fa, fb, d, fd);
    if ((0 == --count) || (fa == 0) || tol.met(a, b)) {
      break;
    }
    prof = (std::fabs(fa - fb) < min_diff) ||
           (std::fabs(fa - fd) < min_diff) ||
           (std::fabs(fa - fe) < min_diff) ||
           (std::fabs(fb - fd) < min_diff) ||
           (std::fabs(fb - fe) < min_diff) ||
           (std::fabs(fd - fe) < min_diff);
    if (prof) {
      c = quadratic_interpolate(a, b, d, fa, fb, fd, 3);
    } else {
      c = cubic_interpolate(a, b, d, e, fa, fb, fd, fe);
    }
    bracket(f, a, b, c, fa, fb, d, fd);
    if ((0 == --count) || (fa == 0) || tol.met(a, b)) {
      break;
    }
    // A double-length secant step from the better end.
    if (std::fabs(fa) < std::fabs(fb)) {
      u = a;
      fu = fa;
    } else {
      u = b;
      fu = fb;
    }
    c = u - 2 * (fu / (fb - fa)) * (b - a);
    if (std::fabs(c - u) > (b - a) / 2) {
      c = a + (b - a) / 2;
    }
    e = d;
    fe = fd;
    bracket(f, a, b, c, fa, fb, d, fd);
    if ((0 == --count) || (fa == 0) || tol.met(a, b)) {
      break;
    }
    // A bisection, where these three steps did not halve the bracket.
    if ((b - a) < mu * (b0 - a0)) {
      continue;
    }
    e = d;
    fe = fd;
    bracket(f, a, b, double(a + (b - a) / 2), fa, fb, d, fd);
    --count;
  }

  max_iter -= count;
  if (fa == 0) {
    b = a;
  } else if (fb == 0) {
    a = b;
  }
  return std::make_pair(a, b);
}

}  // namespace internal

// The bracket TOMS748 ends on, searching [a, b] from f(a) = fa and f(b) = fb,
// which the caller has already evaluated: fa and fb must differ in sign, or
// one be zero, and a < b. f(first) and f(second) differ in sign, or
// first == second where f was exactly zero there. Refuses a bracket that
// misses `tol` after tol.max_iterations evaluations (util::RootFindError), and
// ends that do not bracket a root (util::BracketError).
template <class F>
std::pair<double, double> toms748_bracket(F&& f, double a, double b, double fa,
                                          double fb, const tolerance& tol) {
  std::uintmax_t used = tol.max_iterations;
  const std::pair<double, double> r =
      internal::toms748_solve(f, a, b, fa, fb, tol, used);
  internal::check_converged("toms748", r, tol, used);
  return r;
}

// The same search, returning the final bracket's midpoint.
template <class F>
double toms748(F&& f, double a, double b, double fa, double fb,
               const tolerance& tol) {
  std::uintmax_t used = tol.max_iterations;
  const std::pair<double, double> r =
      internal::toms748_solve(f, a, b, fa, fb, tol, used);
  internal::check_converged("toms748", r, tol, used);
  return (r.first + r.second) / 2.0;
}

// The same, evaluating f at both ends itself. Those two evaluations count
// against tol.max_iterations, so a cap of 2 or less evaluates nothing, and
// returns the midpoint of [a, b] only where [a, b] already meets the
// tolerance. Which end is evaluated first is unspecified, as in Boost: a target
// whose value depends on the order it is called in should take the form above.
template <class F>
double toms748(F&& f, double a, double b, const tolerance& tol) {
  if (a >= b) [[unlikely]] {
    internal::refuse_bracket("toms748", "the ends are out of order", a, b,
                             std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::quiet_NaN());
  }
  std::uintmax_t used = tol.max_iterations;
  std::pair<double, double> r(a, b);
  if (used > 2) {
    used -= 2;
    r = internal::toms748_solve(f, a, b, f(a), f(b), tol, used);
    used += 2;
  }
  internal::check_converged("toms748", r, tol, used);
  return (r.first + r.second) / 2.0;
}

// Bisection of [a, b], evaluating f at both ends itself, which must differ in
// sign (or one be zero) with a < b; returns the final bracket's midpoint. The
// two end evaluations and the first comparison count three against
// tol.max_iterations, as in Boost. For a target with no smoothness to exploit.
// Refuses as toms748 does.
template <class F>
double bisect(F&& f, double a, double b, const tolerance& tol) {
  double min = a;
  double max = b;
  double fmin = f(min);
  double fmax = f(max);
  if (fmin == 0) {
    return min;
  }
  if (fmax == 0) {
    return max;
  }
  if (min >= max) [[unlikely]] {
    internal::refuse_bracket("bisect", "the ends are out of order", a, b,
                             fmin, fmax);
  }
  // A NaN at either end passes, as in Boost.
  if (fmin * fmax >= 0) [[unlikely]] {
    internal::refuse_bracket("bisect", "the ends do not bracket a root", a, b,
                             fmin, fmax);
  }

  std::uintmax_t count = tol.max_iterations;
  if (count < 3) {
    count = 0;
  } else {
    count -= 3;
  }

  while (count && (0 == tol.met(min, max))) {
    double mid = (min + max) / 2;
    double fmid = f(mid);
    if ((mid == max) || (mid == min)) {
      break;
    }
    if (fmid == 0) {
      min = max = mid;
      break;
    } else if (internal::sign(fmid) * internal::sign(fmin) < 0) {
      max = mid;
    } else {
      min = mid;
      fmin = fmid;
    }
    --count;
  }

  const std::pair<double, double> r(min, max);
  internal::check_converged("bisect", r, tol, tol.max_iterations - count);
  return (min + max) / 2;
}

}  // namespace root
}  // namespace odelia

// The functions in `root::internal`, and the searches in toms748 and bisect,
// are adapted from Boost.Math, which carries this notice:
//
//   (C) Copyright John Maddock 2006.
//   (C) Copyright Matt Borland 2024.
//
//   Boost Software License - Version 1.0 - August 17th, 2003
//
//   Permission is hereby granted, free of charge, to any person or
//   organization obtaining a copy of the software and accompanying
//   documentation covered by this license (the "Software") to use,
//   reproduce, display, distribute, execute, and transmit the Software, and
//   to prepare derivative works of the Software, and to permit third-parties
//   to whom the Software is furnished to do so, all subject to the following:
//
//   The copyright notices in the Software and this entire statement,
//   including the above license grant, this restriction and the following
//   disclaimer, must be included in all copies of the Software, in whole or
//   in part, and all derivative works of the Software, unless such copies or
//   derivative works are solely in the form of machine-executable object
//   code generated by a source language processor.
//
//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND
//   NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE
//   DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY,
//   WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
//   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//   SOFTWARE.

#endif
