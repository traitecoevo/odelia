// -*-c++-*-
#ifndef ODELIA_INTERPOLATOR_HPP
#define ODELIA_INTERPOLATOR_HPP

#include <vector>
#include <list>
#include <cmath>
#include <limits>
#include <XAD/XAD.hpp>
#include <odelia/spline.hpp>
#include <odelia/ode_util.hpp>

namespace odelia {
namespace interpolator {

// One interpolator with two build paths: `construct` adaptively refines a node set
// and records it; `init` builds on a node set it is given. A recorded interpolator is
// rebuilt with `init` on those nodes, never re-refined.
//
// Templated on the scalar S of the knot VALUES; knot positions stay double. S = double
// is the production type; an AD active S makes the interpolated value differentiable
// w.r.t. the knot values, delegating to basic_spline<S>.
//
// The basic_interpolator name and the Interpolator alias below are both kept so plant
// compiles unchanged; collapsing them to one Interpolator<S> is done together with
// plant's binding.
template <typename S>
class basic_interpolator {
public:
  // Adaptively refine `target` over [a, b] to tolerance, then build on the chosen
  // nodes. A node is accepted when its midpoint's absolute OR relative error is under
  // tolerance; refinement halves the spacing until every interval passes or the depth
  // cap bites. Decisions are taken in `double` (xad::value), so placement never
  // depends on an active tape -- call this only on the adaptive pass; a replay uses
  // `init` on the recorded nodes.
  template <typename Function>
  void construct(Function target, double a, double b,
                 double atol = 1e-6, double rtol = 1e-6,
                 std::size_t nbase = 17, std::size_t max_depth = 16) {
    if (a >= b) {
      util::stop("Interpolator::construct: impossible bounds (a >= b)");
    }
    if (!util::is_finite(a) || !util::is_finite(b)) {
      util::stop("Interpolator::construct: infinite bounds");
    }
    if (nbase < 2) {
      util::stop("Interpolator::construct: need at least 2 base points");
    }

    // Lists so points can be inserted in the middle of a span during refinement.
    std::list<double> xs;
    std::list<S>      ys;
    std::list<bool>   refine_here;   // is the interval ending at this node still open?

    double dx = (b - a) / static_cast<double>(nbase - 1);
    const double dxmin = dx / std::pow(2.0, static_cast<double>(max_depth));
    for (std::size_t i = 0; i < nbase; ++i) {
      const double xi = a + dx * static_cast<double>(i);
      xs.push_back(xi);
      ys.push_back(target(xi));
      refine_here.push_back(i > 0);
    }
    rebuild(xs, ys);

    auto within_tol = [&](S y_true, S y_pred) {
      const double t = xad::value(y_true), p = xad::value(y_pred);
      return std::fabs(t - p) < atol || std::fabs(1.0 - p / t) < rtol;
    };

    bool open = true;
    while (open) {
      dx /= 2.0;
      if (dx < dxmin) {
        util::stop("Interpolator::construct: refined as far as max_depth allows");
      }
      open = false;
      auto xi = xs.begin();
      auto yi = ys.begin();
      auto zi = refine_here.begin();
      for (; xi != xs.end(); ++xi, ++yi, ++zi) {
        if (*zi) {
          const double x_mid = *xi - dx;
          const S      y_mid = target(x_mid);
          const S      p_mid = eval(x_mid);
          xs.insert(xi, x_mid);
          ys.insert(yi, y_mid);
          const bool still_open = !within_tol(y_mid, p_mid);
          *zi = still_open;                 // the interval [x_mid, *xi]
          refine_here.insert(zi, still_open); // the interval ending at x_mid
          open = open || still_open;
        }
      }
      rebuild(xs, ys);
    }
  }

  // Build an interpolator out of the vectors 'x' and 'y'.
  void init(const std::vector<double> &x_,
            const std::vector<S> &y_) {
    util::check_length(y_.size(), x_.size());
    if (x_.size() < 3)
    {
      util::stop("insufficient number of points");
    }
    x = x_;
    y = y_;
    initialise();
  }

  // Compute the interpolated function from the points contained in 'x' and 'y'.
  void initialise() {
    // https://stackoverflow.com/questions/17769114/stdis-sorted-and-strictly-less-comparison
    if (not std::is_sorted(x.begin(), x.end(), std::less_equal<double>()))
    {
      util::stop("spline control points must be unique and in ascending order");
    }
    if (x.size() > 0)
    {
      spline.set_points(x, y);
      active = true;
    }
  }

  // Support for adding points in turn (assumes monotonic increasing in
  // 'x', unchecked).
  void add_point(double xi, S yi) {
    x.push_back(xi);
    y.push_back(yi);
  }

  // adds point in sorted position (slower than above)
  void add_point_sorted(double xi, S yi) {
    auto x_upper = std::upper_bound(x.begin(), x.end(), xi); // find smallest number larger than xi
    x.insert(x_upper, xi);                                   // add xi below that number
    auto y_upper = std::upper_bound(y.begin(), y.end(), yi);
    y.insert(y_upper, yi);
  }

  // Remove all the contents, being ready to be refilled.
  void clear() {
    x.clear();
    y.clear();
    active = false;
  }

  // Compute the value of the interpolated function at point `x=u`
  S eval(double u) const {
    check_active();
    if (not extrapolate and (u < min() or u > max()))
    {
      util::stop("Extrapolation disabled and evaluation point outside of interpolated domain.");
    }
    return spline(u);
  }

  // faster version of above
  S operator()(double u) const {
    return spline(u);
  }

  // Analytic first derivative dy/du at u (exact derivative of the interpolating
  // polynomial; see Spline::deriv). Useful for exact/smooth gradients.
  S deriv(double u) const {
    check_active();
    return spline.deriv(u);
  }

  // Return the number of (x,y) pairs contained in the Interpolator.
  size_t size() const {
    return x.size();
  }

  // These are chosen so that if a Interpolator is empty, functions
  // looking to see if they will fall outside of the covered range will
  // always find they do.  This is the same principle as R's
  // range(numeric(0)) -> c(Inf, -Inf)
  double min() const {
    return size() > 0 ? x.front() : std::numeric_limits<double>::infinity();
  }

  double max() const {
    return size() > 0 ? x.back() : -std::numeric_limits<double>::infinity();
  }

  void set_extrapolate(bool e) {
    extrapolate = e;
  }

  std::vector<double> get_x() const {
    return x;
  }

  std::vector<S> get_y() const {
    return y;
  }

  // Compute the value of the interpolated function at a vector of
  // points `x=u`, returning a vector of the same length.
  // change to const& vec?
  std::vector<S> r_eval(std::vector<double> u) const {
    check_active();
    auto ret = std::vector<S>();
    ret.reserve(u.size()); // fast to do this once rather than multiple times with push_back
    for (auto const &x : u)
    {
      ret.push_back(eval(x));
    }
    return ret;
  }

private:
  void check_active() const {
    if (!active)
    {
      util::stop("Interpolator not initialised -- cannot evaluate");
    }
  }

  // Rebuild the spline from the working lists during construct().
  void rebuild(const std::list<double>& xs, const std::list<S>& ys) {
    init(std::vector<double>(xs.begin(), xs.end()),
         std::vector<S>(ys.begin(), ys.end()));
  }

  std::vector<double> x;
  std::vector<S> y;
  spline::basic_spline<S> spline;
  bool active = false;
  bool extrapolate = true;
};

// Default interpolator (knot values in double): the production type bound by
// plant's RcppR6 as `odelia::interpolator::Interpolator`, used by ResourceSpline,
// the leaf model, etc.
using Interpolator = basic_interpolator<double>;

}
}

#endif
