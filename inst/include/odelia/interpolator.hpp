// -*-c++-*-
#ifndef ODELIA_INTERPOLATOR_HPP
#define ODELIA_INTERPOLATOR_HPP

#include <vector>
#include <limits>
#include <odelia/spline.hpp>
#include <odelia/ode_util.hpp>

namespace odelia {
namespace interpolator {

// The front end every consumer names: hermite_spline<S> from spline.hpp with
// the calls the family has always made on an interpolator -- init(x, y) from
// values alone, eval with an out-of-domain refusal, deriv, min/max, r_eval --
// added on top. init(x, y) reads the natural cubic spline, as before 0.5.0. Templated on the scalar S of the knot VALUES (knot positions
// stay double); `Interpolator` below pins S = double. A caller that has slopes
// to give uses init(x, y, m), or any of the base class's own calls.
template <typename S>
class hermite_interpolator : public hermite_spline<S> {
  using base = hermite_spline<S>;
public:
  using base::init;   // init(x, y, m): a slope at each knot, supplied

  // From values alone: the natural cubic spline's own knot slopes, so the
  // interpolant is the natural cubic spline odelia has always read here.
  void init(const std::vector<double> &x_, const std::vector<S> &y_) {
    util::check_length(y_.size(), x_.size());
    if (x_.size() < 3)
    {
      util::stop("insufficient number of points");
    }
    check_sorted(x_);
    base::init(x_, y_, natural_slopes<S>(x_, y_));
  }

  // Support for adding points in turn, then initialise(). Assumes increasing x.
  void add_point(double xi, S yi) {
    pending_x.push_back(xi);
    pending_y.push_back(yi);
  }
  void initialise() {
    init(pending_x, pending_y);
    pending_x.clear();
    pending_y.clear();
  }

  // Remove all the contents, being ready to be refilled.
  void clear() {
    pending_x.clear();
    pending_y.clear();
    *static_cast<base*>(this) = base();
  }

  // Compute the value of the interpolated function at point `x=u`
  S eval(double u) const {
    check_active();
    // ⚠️ Do NOT "tidy" this into `not (u >= min() and u <= max())`. Every
    // comparison against NaN is false, so as written a non-finite `u` falls
    // *through* to the spline and comes back non-finite -- which callers rely on
    // (traitecoevo/plant#576 documents a `profit_psi_stem_TF(NA, .) -> NA`
    // contract built on it). Negating an in-range test turns that into a throw:
    // it reads as a tightening and is a behaviour change.
    if (not extrapolate and (u < min() or u > max()))
    {
      const bool below = u < min();
      // The point, how far out it fell, and the domain. Reporting none of the
      // three used to make an out-of-domain failure a bisect rather than a read:
      // localising plant#576 meant instrumenting four call sites by hand to
      // discover which spline was being asked and at what value, and the answer
      // (the LOWER end, not past the far end as everyone assumed) inverted the
      // fix. The caller's own identity is the one thing this layer cannot know --
      // consumers that build several splines should catch and say which.
      util::stop(std::string("Extrapolation disabled and evaluation point "
                             "outside of interpolated domain: u = ") +
                 util::format_double(u) + " lies " +
                 util::format_double(below ? min() - u : u - max()) +
                 " beyond the " + (below ? "lower" : "upper") + " end of [" +
                 util::format_double(min()) + ", " +
                 util::format_double(max()) + "].");
    }
    return base::eval(u);
  }

  // eval() without its checks: no domain refusal and no initialisation check, so
  // reading an empty interpolator is undefined. For hot loops whose caller has
  // already bounded u, as plant's light field does per quadrature point.
  S operator()(double u) const {
    return base::eval_unchecked(u);
  }

  // Analytic first derivative dy/du at u (exact derivative of the interpolating
  // polynomial). Useful for exact/smooth gradients.
  S deriv(double u) const {
    check_active();
    return base::slope(u);
  }

  // These are chosen so that if a Interpolator is empty, functions
  // looking to see if they will fall outside of the covered range will
  // always find they do.  This is the same principle as R's
  // range(numeric(0)) -> c(Inf, -Inf)
  double min() const {
    return this->size() > 0 ? this->knots().front()
                            : std::numeric_limits<double>::infinity();
  }

  double max() const {
    return this->size() > 0 ? this->knots().back()
                            : -std::numeric_limits<double>::infinity();
  }

  void set_extrapolate(bool e) {
    extrapolate = e;
  }

  std::vector<double> get_x() const {
    return this->knots();
  }

  std::vector<S> get_y() const {
    return this->values();
  }

  // Compute the value of the interpolated function at a vector of
  // points `x=u`, returning a vector of the same length.
  std::vector<S> r_eval(std::vector<double> u) const {
    check_active();
    auto ret = std::vector<S>();
    ret.reserve(u.size());
    for (auto const &x : u)
    {
      ret.push_back(eval(x));
    }
    return ret;
  }

private:
  static void check_sorted(const std::vector<double> &x) {
    // https://stackoverflow.com/questions/17769114/stdis-sorted-and-strictly-less-comparison
    if (not std::is_sorted(x.begin(), x.end(), std::less_equal<double>()))
    {
      util::stop("spline control points must be unique and in ascending order");
    }
  }

  void check_active() const {
    if (this->size() == 0)
    {
      util::stop("Interpolator not initialised -- cannot evaluate");
    }
  }

  std::vector<double> pending_x;
  std::vector<S> pending_y;
  bool extrapolate = true;
};

// Default interpolator (knot values in double): the production type used by
// the drivers, plant's ResourceSpline, the leaf model, etc.
using Interpolator = hermite_interpolator<double>;

}
}

#endif
