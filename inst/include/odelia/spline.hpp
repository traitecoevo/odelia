// -*-c++-*-
#ifndef ODELIA_SPLINE_HPP
#define ODELIA_SPLINE_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <list>
#include <type_traits>
#include <limits>
#include <vector>
#include <odelia/ode_util.hpp>

// Interpolation: one interpolant, and the two rules that produce its inputs.
//
// A piecewise cubic is determined by a value and a slope at each knot, so the
// interpolant takes exactly those and holds no rule for choosing any of them.
// Where the knot data comes from belongs to the caller, and there are three
// sources: a closed form supplies both halves, a reduction supplies both halves
// from one expression, or only values exist and `natural_slopes` (the default)
// or `monotone_slopes` chooses the
// rest. Where the knot POSITIONS come from belongs to the caller too: a fixed
// lattice, or `refine` for a target whose features are not known in advance.
//
// There is no fit that chooses slopes globally. One was here, solving a
// tridiagonal system for C2, and it converged at h^2 on the curve this library
// tabulates against h^3.7 for the same knots read as a Hermite -- the global
// solve spreads a local defect -- while making every knot influence every span,
// which turns an O(1) adjoint into an O(K) one. Nothing read a second derivative
// from it.
namespace odelia {
namespace interpolator {

// A piecewise-polynomial interpolant built from a value AND a slope at each knot,
// and a curvature as well where the source has one.
//
// The value and the slope come from one polynomial, so a caller that needs both
// gets a consistent pair: slope(u) is the exact derivative of what eval(u) returns.
//
// Order is how many coefficients a span holds, and it is set by which set_data a
// caller reaches for: two channels build a cubic, three build a quintic, and there
// is one polynomial written here for both. A cubic converges as h^4 and a quintic
// as h^6, so a source with a closed form for the second derivative reaches a given
// error on far fewer knots -- and a quantity that is not twice differentiable
// simply has no third channel to supply.
//
// ⚠️ IT IS NOT A SETTING TO BE COLLAPSED, and the reason is the tape rather than the
// clock. Reading every span as a quintic costs read time, which is arguable --
// but a cubic's top two coefficients would then be ACTIVE zeros holding
// tape slots, so each read would record two operations whose adjoints are
// structurally zero. That doubles nothing and adds half again to the sparsest,
// hottest recorded read in the model.
//
// Each span reads only its own two knots, so moving one knot changes the
// interpolant only in the two spans that touch it.
//
// Knot positions are double; values and slopes carry the working scalar S. The two
// halves of a build are separate: set_nodes lays out the spans from the positions,
// and set_data fills the coefficients. A caller whose positions are fixed for a run
// calls set_nodes once and set_data per stage.
//
// eval takes either a double position or an active one; at an active position the
// value is read at its passive part and the query's own derivative is carried on
// through the slope, so d(value)/d(u) is recorded. slope and value_and_slope take a
// double, because a query's derivative reaches a value and not a slope.
template <typename S, int Order = 3>
class hermite_spline {
  static_assert(Order == 3 || Order == 5,
                "a span is determined by a value and a slope at each end, or by "
                "those plus a curvature; there is no other knot data any source "
                "in this library supplies exactly.");

public:
  // Knot positions, strictly ascending; at least two are needed (one span).
  // Discards any data already set.
  void set_nodes(const std::vector<double>& x_) {
    if (x_.size() < 2) util::stop("hermite_spline: need at least 2 knots");
    for (std::size_t i = 1; i < x_.size(); ++i) {
      if (!(x_[i] > x_[i - 1]))
        util::stop("hermite_spline: knots must be strictly ascending");
    }
    x = x_;
    const std::size_t ns = x.size() - 1;
    spans.assign(ns, Span());
    for (std::size_t k = 0; k < ns; ++k) {
      spans[k].x0 = x[k];
      spans[k].inv_h = 1.0 / (x[k + 1] - x[k]);
    }
    // An equally spaced knot set indexes by arithmetic instead of a search.
    uniform = false;
    if (ns > 1) {
      const double h0 = x[1] - x[0];
      const double tol = 1e-12 * (x.back() - x.front());
      uniform = true;
      for (std::size_t k = 1; k < ns; ++k) {
        if (std::abs((x[k + 1] - x[k]) - h0) > tol) { uniform = false; break; }
      }
      if (uniform) inv_h0 = 1.0 / h0;
    }
    inv_mean_h = static_cast<double>(ns) / (x.back() - x.front());
    initialised = false;
  }

  // Values and dy/dx at the nodes already set, one entry each per node.
  void set_data(const std::vector<S>& y_, const std::vector<S>& dydx_) {
    static_assert(Order == 3,
                  "a quintic span has a curvature coefficient, so it needs the "
                  "second derivative at every knot; leaving it out would make the "
                  "read quintic in storage and cubic in fact.");
    begin_data(y_, dydx_);
    for (std::size_t k = 0; k < spans.size(); ++k) {
      const double h = x[k + 1] - x[k];
      const S a = y[k], b = y[k + 1];
      const S sa = m[k] * h, sb = m[k + 1] * h;
      Span& s = spans[k];
      s.y0 = a;
      s.c[0] = sa;
      s.c[1] = 3.0 * (b - a) - 2.0 * sa - sb;
      s.c[2] = 2.0 * (a - b) + sa + sb;
    }
    initialised = true;
  }

  // As above, with d2y/dx2 as well: the unique quintic through both knots' value,
  // slope and curvature.
  void set_data(const std::vector<S>& y_, const std::vector<S>& dydx_,
                const std::vector<S>& d2ydx2_) {
    static_assert(Order == 5,
                  "a cubic span is determined by the value and the slope, so a "
                  "curvature has nowhere to go and would be silently dropped.");
    begin_data(y_, dydx_);
    util::check_length(d2ydx2_.size(), x.size());
    // The curvatures are read into the spans and not kept. values() and slopes()
    // are held because a caller hands them on and a slope recovered from a span is
    // not bit-identical; nothing hands a curvature on, and a member here would sit
    // on every cubic interpolant in the library.
    for (std::size_t j = 0; j < spans.size(); ++j) {
      const double h = x[j + 1] - x[j];
      const S a = y[j], b = y[j + 1];
      const S a1 = m[j] * h, b1 = m[j + 1] * h;
      const S a2 = 0.5 * d2ydx2_[j] * h * h, b2 = d2ydx2_[j + 1] * h * h;
      // What the two ends leave for the top three powers to match.
      const S d = b - a - a1 - a2;
      const S e = b1 - a1 - 2.0 * a2;
      const S f = b2 - 2.0 * a2;
      Span& s = spans[j];
      s.y0 = a;
      s.c[0] = a1;
      s.c[1] = a2;
      s.c[2] = 10.0 * d - 4.0 * e + 0.5 * f;
      s.c[3] = -15.0 * d + 7.0 * e - f;
      s.c[4] = 6.0 * d - 3.0 * e + 0.5 * f;
    }
    initialised = true;
  }

  // Nodes and data in one call, for a caller that rebuilds both together.
  void init(const std::vector<double>& x_, const std::vector<S>& y_,
            const std::vector<S>& dydx_) {
    set_nodes(x_);
    set_data(y_, dydx_);
  }

  void init(const std::vector<double>& x_, const std::vector<S>& y_,
            const std::vector<S>& dydx_, const std::vector<S>& d2ydx2_) {
    set_nodes(x_);
    set_data(y_, dydx_, d2ydx2_);
  }

  void clear() {
    x.clear(); y.clear(); m.clear(); spans.clear();
    inv_h0 = 0.0;
    inv_mean_h = 0.0;
    uniform = false;
    initialised = false;
  }

  bool is_initialised() const { return initialised; }
  std::size_t size() const { return x.size(); }
  // ⚠️ NOT GUARDED, and a consumer that may read a domain before filling one owes
  // its own check: these are `x.front()` and `x.back()` on an empty vector.
  // `eval` checks, because it can afford to; these two sit inside a caller's own
  // bound test, which is where the hot path reads them.
  double min() const { return x.front(); }
  double max() const { return x.back(); }
  const std::vector<double>& knots() const { return x; }

  // The data the last set_data was given, unchanged. A caller that needs to hand
  // the interpolant's contents on -- to pack them into a state vector, or to
  // supply them to a pass that does not build them -- reads them here rather than
  // keeping its own copy: a slope recovered from a span is not bit-identical,
  // where these are what was supplied.
  const std::vector<S>& values() const { return y; }
  const std::vector<S>& slopes() const { return m; }

  // Value at u. Outside the knot range the end span's line is extended (value and
  // slope of the nearest end), which keeps the read C1 across the boundary instead
  // of letting a cubic run away. A caller that must refuse an out-of-range read
  // compares against min() and max() and says which curve it was asking about,
  // which this class cannot know.
  template <typename U>
  S eval(const U& u) const {
    check_initialised();
    const double up = util::to_passive(u);
    if constexpr (std::is_same_v<U, double>) {
      return value_at(up);
    } else {
      // One span load for both halves. The with_query_derivative needs the slope too, and taking
      // them separately resolves the same position twice on the read this class
      // exists for.
      S value, dydu;
      value_and_slope(up, value, dydu);
      return with_query_derivative(value, dydu, u, up);
    }
  }

  template <typename U>
  S operator()(const U& u) const { return eval(u); }

protected:
  // eval at a double position without the initialisation check, for a front end
  // whose own unchecked read is documented as such.
  S eval_unchecked(double u) const { return value_at(u); }

public:

  // dy/du at u -- the exact derivative of the polynomial eval() uses, as a value.
  //
  // The position is read passively and an active query is refused rather than
  // answered with a silent zero. Carrying the query's derivative here would need
  // the span's second derivative, and that number is a derivative of the fit: the
  // read is C1 and not C2, so a curvature taken from it describes the interpolant
  // rather than what was interpolated.
  template <typename U>
  S slope(const U& u) const {
    static_assert(std::is_same_v<U, double>,
                  "hermite_spline::slope reads the position at its value, so "
                  "an active query's derivative has nowhere to go and would come "
                  "back as exactly zero. eval() is the reader that takes one.");
    check_initialised();
    return slope_at(u);
  }

  // Both from one knot lookup and one span load, so a caller wanting the pair at
  // many positions pays one lookup each rather than two.
  //
  // Passive query only, for the reason slope() is: the pair's second half has no
  // route for the query's derivative, so carrying it onto the first half alone
  // would answer half the query and leave the other half reading exactly zero.
  template <typename U>
  void value_and_slope(const U& u, S& value, S& dydu) const {
    static_assert(std::is_same_v<U, double>,
                  "hermite_spline::value_and_slope reads the position at its "
                  "value, so an active query's derivative would reach the value and "
                  "not the slope. eval() is the reader that takes one.");
    check_initialised();
    const double up = u;
    if (up <= x.front()) {
      value = y.front() + m.front() * (up - x.front());
      dydu = m.front();
    } else if (up >= x.back()) {
      value = y.back() + m.back() * (up - x.back());
      dydu = m.back();
    } else {
      const Span& s = spans[span_of(up)];
      const double t = s.local(up);
      value = s.value(t);
      dydu = s.slope(t);
    }
  }

private:
  // One span's whole polynomial, contiguous: a query touches a single cache line
  // rather than one per coefficient array. The cubic and its derivative are written
  // here, beside the coefficients they read, so the three readers share one spelling
  // of each instead of restating it and needing a test that they agree.
  struct Span {
    double x0 = 0.0, inv_h = 0.0;
    S y0{};
    S c[Order]{};  // the coefficient of t^(i+1), t the span-local coordinate
    double local(double u) const { return (u - x0) * inv_h; }
    // Horner over the coefficients the order has, once for both orders. Order is a
    // constant so both loops unroll, and the factors stay double: at an active
    // scalar a converted int would register a variable and record an operation
    // per term.
    S value(double t) const {
      S r = c[Order - 1];
      for (int i = Order - 2; i >= 0; --i) r = c[i] + t * r;
      return y0 + t * r;
    }
    S slope(double t) const {
      S r = static_cast<double>(Order) * c[Order - 1];
      for (int i = Order - 2; i >= 1; --i)
        r = static_cast<double>(i + 1) * c[i] + t * r;
      // c[0]'s factor is one, and writing it as one is a multiply the compiler
      // keeps: it measured 8% of the read.
      return (c[0] + t * r) * inv_h;
    }
  };

  // The two length checks and the two vectors every order stores.
  void begin_data(const std::vector<S>& y_, const std::vector<S>& dydx_) {
    if (spans.empty()) util::stop("hermite_spline: no knots set");
    util::check_length(y_.size(), x.size());
    util::check_length(dydx_.size(), x.size());
    y = y_;
    m = dydx_;
  }

  // The query's derivative, materialised while its operands are alive. A deduced
  // return type here would hand back an XAD expression template referencing the
  // temporaries of this return statement, which die on return.
  template <typename U>
  static S with_query_derivative(const S& value, const S& dydu, const U& u, double up) {
    static_assert(std::is_constructible_v<S, U>,
                  "hermite_spline: reading at an active position needs the "
                  "knot values on the same scalar, so the derivative of the query "
                  "has somewhere to go -- an active position with S = double would "
                  "silently drop it.");
    return value + dydu * (u - up);
  }

  S value_at(double u) const {
    if (u <= x.front()) return y.front() + m.front() * (u - x.front());
    if (u >= x.back())  return y.back()  + m.back()  * (u - x.back());
    const Span& s = spans[span_of(u)];
    return s.value(s.local(u));
  }

  S slope_at(double u) const {
    if (u <= x.front()) return m.front();
    if (u >= x.back())  return m.back();
    const Span& s = spans[span_of(u)];
    return s.slope(s.local(u));
  }

  std::size_t span_of(double u) const {
    const std::size_t ns = spans.size();
    // A non-finite query reaches here only as NaN -- the infinities are answered by
    // the end-extension branches above -- and casting NaN to an index is undefined.
    // Answer from a span that exists and let the arithmetic carry the NaN out, so
    // a non-finite position reads back non-finite rather than throwing.
    if (!util::is_finite(u)) return 0;
    if (uniform) {
      const std::size_t k = static_cast<std::size_t>((u - x.front()) * inv_h0);
      return k < ns ? k : ns - 1;
    }
    // A graded grid: guess from the mean spacing, then step to the span with
    // x[k] <= u < x[k + 1] -- the span a binary search returns, so the read is
    // unchanged. Every caller has answered u <= x.front() and u >= x.back()
    // already, so u is strictly inside.
    //
    // ⚠️ Do not replace this with std::upper_bound, and do not cap the steps with a
    // binary-search fallback. Queries typically arrive in order, so the walk from
    // the guess is short and either replacement is slower; a synthetic benchmark
    // with random queries favours the search and is the wrong test.
    std::size_t k = static_cast<std::size_t>((u - x.front()) * inv_mean_h);
    if (k > ns - 1) k = ns - 1;
    while (k > 0 && x[k] > u) --k;
    while (k + 1 < ns && x[k + 1] <= u) ++k;
    return k;
  }

  void check_initialised() const {
    if (!initialised) util::stop("hermite_spline: not initialised");
  }

  std::vector<double> x;   // knot positions, contiguous for the search
  std::vector<S> y, m;     // knot values and slopes, as supplied
  std::vector<Span> spans;

public:
  // Every value here that carries the working scalar. The knot positions are not
  // among them: they stay double, which is what keeps the grid out of the
  // derivative.
  // Both overloads, for the reason value_with_slope.hpp gives: visit_active dispatches
  // on whether the call compiles, so a const object handed to a non-const walk is
  // passed over in silence rather than refused, and the rewinding forms in
  // implicit_node.hpp take their inputs const.
  template <class F>
  void for_each_active(F&& f) {
    for (S& v : y) { f(v); }
    for (S& v : m) { f(v); }
    for (Span& span : spans) {
      f(span.y0);
      for (S& coefficient : span.c) { f(coefficient); }
    }
  }
  template <class F>
  void for_each_active(F&& f) const {
    for (const S& v : y) { f(v); }
    for (const S& v : m) { f(v); }
    for (const Span& span : spans) {
      f(span.y0);
      for (const S& coefficient : span.c) { f(coefficient); }
    }
  }

private:
  double inv_h0 = 0.0;
  double inv_mean_h = 0.0;  // spans per unit of the knot range, for the lookup's guess
  bool uniform = false;
  bool initialised = false;
};

// The knot slopes of the natural cubic spline through (x, y): the slopes a
// global C2 fit with zero curvature at both ends would give. A cubic Hermite
// read through these slopes IS that spline on every span, so this is the rule
// that reproduces the interpolant odelia had before 0.5.0, and it is the default
// for a caller with values alone. Solved once, at build, by the tridiagonal
// (Thomas) recurrence on the half-curvatures b_i, b_0 = b_{n-1} = 0:
//
//   h_{i-1}/3 b_{i-1} + 2(h_{i-1} + h_i)/3 b_i + h_i/3 b_{i+1}
//       = (y_{i+1} - y_i)/h_i - (y_i - y_{i-1})/h_{i-1}
//
// and the slope at knot i is (y_{i+1} - y_i)/h_i - (2 b_i + b_{i+1}) h_i / 3,
// with the last knot's taken from the final span's own cubic.
//
// Every slope depends on every value, so at an active scalar a read through
// these slopes reaches every knot; a caller that can supply its own slopes
// should, and gets a read that reaches four.
template <class S>
inline std::vector<S> natural_slopes(const std::vector<double>& x,
                                     const std::vector<S>& y) {
  util::check_length(y.size(), x.size());
  const std::size_t n = x.size();
  if (n < 2) util::stop("natural_slopes: need at least 2 knots");
  std::vector<S> m(n);
  if (n == 2) {
    m[0] = m[1] = (y[1] - y[0]) / (x[1] - x[0]);
    return m;
  }
  std::vector<double> h(n - 1);
  for (std::size_t i = 0; i + 1 < n; ++i) h[i] = x[i + 1] - x[i];
  // Forward sweep over the interior rows; c_prime and d_prime are the
  // eliminated super-diagonal and right-hand side.
  std::vector<double> c_prime(n, 0.0);
  std::vector<S> d_prime(n, S(0.0));
  for (std::size_t i = 1; i + 1 < n; ++i) {
    const double lower = h[i - 1] / 3.0;
    const double diag = 2.0 * (h[i - 1] + h[i]) / 3.0;
    const double upper = h[i] / 3.0;
    const S rhs = (y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1];
    const double denom = diag - lower * c_prime[i - 1];
    c_prime[i] = upper / denom;
    d_prime[i] = (rhs - lower * d_prime[i - 1]) / denom;
  }
  std::vector<S> b(n, S(0.0));
  for (std::size_t i = n - 2; i >= 1; --i) {
    b[i] = d_prime[i] - c_prime[i] * b[i + 1];
  }
  for (std::size_t i = 0; i + 1 < n; ++i) {
    m[i] = (y[i + 1] - y[i]) / h[i] - (2.0 * b[i] + b[i + 1]) * h[i] / 3.0;
  }
  const double hl = h[n - 2];
  const S a_last = (b[n - 1] - b[n - 2]) / (3.0 * hl);
  m[n - 1] = 3.0 * a_last * hl * hl + 2.0 * b[n - 2] * hl + m[n - 2];
  return m;
}

// Slopes for knots that arrive with values and nothing else, limited so that the
// fit is monotone on every span -- so a read can never leave the range of the two
// values that bound it.
//
// Double only, and that is the point rather than a limitation: a caller holding
// values without slopes is holding data from outside the model, and data is
// passive. The limit also decides its branch by comparing values, so at an active
// scalar the rule would be differentiable everywhere except where that branch
// switches, which is a derivative no caller here wants.
//
// The interior estimate is the parabola through three values, which on an uneven
// grid is more accurate than the arithmetic mean of the two secants; the limit is
// Fritsch and Carlson, projecting the pair of end slopes back onto the region
// where the span is monotone.
//
// ⚠️ THE REGION IS NOT THE CIRCLE alpha^2 + beta^2 <= 9. That circle sits inside
// it, so testing against it fires on spans that were already monotone and flattens
// them, costing accuracy on smooth data. Both keep an intermittent series inside
// its own values.
template <class S>
inline std::vector<S> monotone_slopes(const std::vector<double>& x,
                                      const std::vector<S>& y) {
  util::check_length(y.size(), x.size());
  if (x.size() < 2) util::stop("monotone_slopes: need at least 2 knots");
  const std::size_t n = x.size();
  std::vector<S> secant(n - 1), m(n);
  for (std::size_t k = 0; k + 1 < n; ++k) {
    secant[k] = (y[k + 1] - y[k]) / (x[k + 1] - x[k]);
  }
  m[0] = secant[0];
  m[n - 1] = secant[n - 2];
  for (std::size_t k = 1; k + 1 < n; ++k) {
    const double h0 = x[k] - x[k - 1], h1 = x[k + 1] - x[k];
    m[k] = (h1 * secant[k - 1] + h0 * secant[k]) / (h0 + h1);
  }
  // A turning point in the data is a turning point in the fit. Without this the
  // estimate at a peak lies between two secants of opposite sign, so it opposes one
  // of them -- and the projection below cannot reach that case, which is how a
  // rainfall series that is nowhere negative reads to -1.86 between two wet days.
  // It costs accuracy where the data really is smooth: on a sine at 100 knots the
  // mean relative difference goes from 1.2e-05 to 4.5e-05. Both are far below the
  // precision of any series a caller supplies, and one of them is a value the model
  // cannot have.
  for (std::size_t k = 1; k + 1 < n; ++k) {
    if (secant[k - 1] * secant[k] <= 0.0) {
      m[k] = 0.0;
    }
  }
  for (std::size_t k = 0; k + 1 < n; ++k) {
    const S s = secant[k];
    if (s == 0.0) {
      // A flat pair pins both its slopes, which is what stops a pulse being
      // smeared back across the dry interval beside it.
      m[k] = 0.0;
      m[k + 1] = 0.0;
      continue;
    }
    const S alpha = m[k] / s, beta = m[k + 1] / s;
    const S a2b3 = 2.0 * alpha + beta - 3.0;
    const S ab23 = alpha + 2.0 * beta - 3.0;
    if (a2b3 > 0.0 && ab23 > 0.0 && alpha * (a2b3 + ab23) < a2b3 * a2b3) {
      using std::sqrt;
      const S tau = 3.0 * s / sqrt(alpha * alpha + beta * beta);
      m[k] = tau * alpha;
      m[k + 1] = tau * beta;
    }
  }
  return m;
}

// A node set and the interpolant's data on it, which is what `refine` returns and
// what `hermite_spline::init` takes.
template <typename S>
struct nodes_and_data {
  std::vector<double> x;
  std::vector<S> y, m;
};

// Choose knot positions for a target whose features are not known in advance:
// halve the spacing until every interval's midpoint is within tolerance of the
// interpolant built on the interval's ends.
//
// The target returns a value AND a slope, which is what the interpolant takes, so
// refinement needs no rule of its own for the second half. Placement is decided
// in double (util::to_passive), so a refined node set never depends on an active
// tape -- and a pass that must run on the nodes a previous pass chose calls
// `hermite_spline::init` on them rather than refining again.
template <typename S, typename Function>
nodes_and_data<S> refine(Function value_and_slope, double a, double b,
                         double tol = 1e-6, std::size_t nbase = 17,
                         std::size_t max_depth = 16) {
  if (!(a < b)) util::stop("interpolator::refine: impossible bounds (a >= b)");
  if (!util::is_finite(a) || !util::is_finite(b)) {
    util::stop("interpolator::refine: infinite bounds");
  }
  if (nbase < 2) util::stop("interpolator::refine: need at least 2 base points");

  // Lists so a node can be inserted in the middle of a span during refinement.
  std::list<double> xs;
  std::list<S> ys, ms;
  std::list<bool> open;   // is the interval ending at this node still open?

  double dx = (b - a) / static_cast<double>(nbase - 1);
  const double dxmin = dx / std::pow(2.0, static_cast<double>(max_depth));
  for (std::size_t i = 0; i < nbase; ++i) {
    const auto vs = value_and_slope(a + dx * static_cast<double>(i));
    xs.push_back(a + dx * static_cast<double>(i));
    ys.push_back(vs.first);
    ms.push_back(vs.second);
    open.push_back(i > 0);
  }

  auto gather = [&]() {
    return nodes_and_data<S>{std::vector<double>(xs.begin(), xs.end()),
                             std::vector<S>(ys.begin(), ys.end()),
                             std::vector<S>(ms.begin(), ms.end())};
  };
  hermite_spline<S> fit;
  { const auto d = gather(); fit.init(d.x, d.y, d.m); }

  bool any_open = true;
  while (any_open) {
    dx /= 2.0;
    if (dx < dxmin) {
      util::stop("interpolator::refine: refined as far as max_depth allows");
    }
    any_open = false;
    auto xi = xs.begin();
    auto yi = ys.begin();
    auto mi = ms.begin();
    auto oi = open.begin();
    for (; xi != xs.end(); ++xi, ++yi, ++mi, ++oi) {
      if (!*oi) continue;
      const double x_mid = *xi - dx;
      const auto vs = value_and_slope(x_mid);
      const double got = util::to_passive(vs.first);
      const double got_fit = util::to_passive(fit.eval(x_mid));
      const bool still_open =
          std::fabs(got - got_fit) > tol * std::max(std::fabs(got), 1.0);
      xs.insert(xi, x_mid);
      ys.insert(yi, vs.first);
      ms.insert(mi, vs.second);
      *oi = still_open;                 // the interval [x_mid, *xi]
      open.insert(oi, still_open);      // the interval ending at x_mid
      any_open = any_open || still_open;
    }
    const auto d = gather();
    fit.init(d.x, d.y, d.m);
  }
  return gather();
}

}
}

#endif
