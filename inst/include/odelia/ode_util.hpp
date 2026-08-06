// -*-c++-*-
#ifndef ODELIA_UTIL_HPP_
#define ODELIA_UTIL_HPP_

// This header, and therefore the whole solver core that includes it, is free of
// R: `stop` throws a std::runtime_error, which Rcpp converts into an ordinary R
// error at the package boundary, and `warning` writes to std::cerr. That is what
// lets consumers of inst/include/ compile and run as plain C++ with no R
// installation at all (traitecoevo/leaf_cpp#11). R belongs in src/ and in the
// interface headers -- solver_interface.hpp, rcpp_interface_helpers.hpp -- not
// here.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <stddef.h> // size_t

namespace odelia {
namespace util {

inline bool is_finite(double x) {
  return std::isfinite(x);
}

// Throws; never returns. The attribute lets the compiler see that callers whose
// error branches end in util::stop() do not fall through.
[[noreturn]] inline void stop(const std::string &msg) {
  throw std::runtime_error(msg);
}

// A state outside the system's valid domain, as distinct from a bug (#55).
//
// The adaptive stepper catches *this type specifically* and treats it as a step
// rejection -- shrink and retry -- rather than letting it end the solve. Anything
// else, util::stop() included, still propagates, and that distinction is
// load-bearing: catching runtime_error broadly would turn a genuine programming
// error into step-shrinking until "Cannot achieve the desired accuracy", a
// diagnostic pointing at the solver instead of at the bug. A system saying "this
// state is impossible" and a system saying "this code is wrong" must not be the
// same signal.
//
// Deriving from std::runtime_error keeps the Rcpp conversion to an ordinary R
// error intact for the cases that do escape -- a fixed-step solve, which has no
// step size to shrink, or an adaptive one already at step_size_min.
struct DomainError : std::runtime_error {
  explicit DomainError(const std::string &msg) : std::runtime_error(msg) {}
};

// As stop(), but for a state the model has no meaning for rather than an error in
// the code. Prefer a message naming the quantity, its value and the constraint:
// this string is what the solver reports if the step cannot be shrunk far enough.
[[noreturn]] inline void stop_domain(const std::string &msg) {
  throw DomainError(msg);
}

// Not an R warning: nothing in the solver core may assume an R session exists.
// Callers that need one should raise it from their own R-facing code. Uses
// fprintf rather than std::cerr to keep <iostream>, and its per-translation-unit
// static initialiser, out of a header that everything includes.
inline void warning(const std::string &msg) {
  std::fprintf(stderr, "odelia: %s\n", msg.c_str());
}

// A double in an error message. std::to_string is fixed-point with six decimals,
// so it renders a flux of 1e-22 as "0.000000" and a domain endpoint of 1e8 with
// eight useless digits -- in both cases erasing the number the reader needed.
// Six significant figures, matching plant's util::format_double so the family
// renders numbers the same way; enough to identify a value, and short enough that
// 6.8918 does not arrive as 6.8917999999999999. Deliberately NOT round-trip
// precision: these strings are for reading, not for reconstructing a double.
inline std::string format_double(double x) {
  char buf[32];
  std::snprintf(buf, sizeof(buf), "%.6g", x);
  return std::string(buf);
}

inline void check_length(size_t received, size_t expected)
{
  if (expected != received)
  {
    stop("Incorrect length input; expected " +
         std::to_string(expected) + ", received " +
         std::to_string(received));
  }
}

// Use this to be explicit when a potentially unsafe floating point
// equality test is being made.
inline
bool identical(double a, double b) {
  std::uint64_t ua = 0;
  std::uint64_t ub = 0;
  std::memcpy(&ua, &a, sizeof(double));
  std::memcpy(&ub, &b, sizeof(double));
  return ua == ub;
}

// http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template<class T>
bool almost_equal(T x, T y, int ulp) {
  // the machine epsilon has to be scaled to the magnitude of the larger value
  // and multiplied by the desired precision in ULPs (units in the last place)
  return std::abs(x - y) <=   std::numeric_limits<T>::epsilon()
    * std::max(std::abs(x), std::abs(y))
    * ulp;
}

// Based on C++11's is_sorted
template <class ForwardIterator>
bool is_sorted(ForwardIterator first, ForwardIterator last) {
  if (first == last)
    return true;

  ForwardIterator next = first;
  while (++next != last) {
    if (*next < *first)
      return false;
    ++first;
  }
  return true;
}

template <class ForwardIterator>
bool is_decreasing(ForwardIterator first, ForwardIterator last) {
  if (first == last)
    return true;

  ForwardIterator next = first;
  while (++next != last) {
    if (*next > *first)
      return false;
    ++first;
  }
  return true;
}

template<typename T>
std::string to_string(T x) {
  return std::to_string(x);
}

}
}

#endif
