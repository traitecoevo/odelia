// -*-c++-*-
#ifndef ODELIA_INTERPOLATOR_HPP
#define ODELIA_INTERPOLATOR_HPP

#include <vector>
#include <odelia/spline.hpp>
#include <odelia/ode_util.hpp>

namespace odelia {
namespace interpolator {

class Interpolator {
public:
  // Build an interpolator out of the vectors 'x' and 'y'.
  void init(const std::vector<double> &x_,
            const std::vector<double> &y_) {
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
  void add_point(double xi, double yi) {
    x.push_back(xi);
    y.push_back(yi);
  }

  // adds point in sorted position (slower than above)
  void add_point_sorted(double xi, double yi) {
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
  double eval(double u) const {
    check_active();
    if (not extrapolate and (u < min() or u > max()))
    {
      util::stop("Extrapolation disabled and evaluation point outside of interpolated domain.");
    }
    return spline(u);
  }

  // faster version of above
  double operator()(double u) const {
    return spline(u);
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
    return size() > 0 ? x.front() : R_PosInf;
  }

  double max() const {
    return size() > 0 ? x.back() : R_NegInf;
  }

  void set_extrapolate(bool e) {
    extrapolate = e;
  }

  std::vector<double> get_x() const {
    return x;
  }

  std::vector<double> get_y() const {
    return y;
  }

  // Compute the value of the interpolated function at a vector of
  // points `x=u`, returning a vector of the same length.
  // change to const& vec?
  std::vector<double> r_eval(std::vector<double> u) const {
    check_active();
    auto ret = std::vector<double>();
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

  std::vector<double> x, y;
  spline::Spline spline;
  bool active = false;
  bool extrapolate = true;
};

}
}

#endif
