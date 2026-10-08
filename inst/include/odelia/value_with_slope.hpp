// -*-c++-*-
#ifndef ODELIA_VALUE_WITH_SLOPE_HPP_
#define ODELIA_VALUE_WITH_SLOPE_HPP_

namespace odelia {

// A quantity and its derivative along whatever the caller is differentiating.
//
// The pair rather than two scalars, because the two are only meaningful together:
// a consumer handed a value and a slope from separate places can pair them across
// different points, different orders, or different independent variables, and all
// three compile. Wherever a value is carried with its own derivative -- a knot on
// an interpolated field, a coordinate at an operating point, a frame of a
// reduction -- this is the type.
//
// WHAT THE SLOPE IS WITH RESPECT TO is the caller's, and this type does not name
// it: one caller's is a height, another's a potential. Naming the variable here
// would make one of them wrong.
//
// It lives here because of for_each_active, not because it is shared. visit_active
// does not open an aggregate of two scalars, so a pair that does not say what it
// holds loses both members from the walk in silence (see visit_active). That
// obligation is this library's, so a model defining the pair itself is a model
// carrying this library's problem.
//
// No includes of its own, so a header on any path can take it without taking
// anything else with it.
template <typename T>
struct value_with_slope {
  T value;
  T slope;

  // Both, because visit_active refuses a type that declares the non-const walk
  // without the const one.
  template <class F>
  void for_each_active(F&& f) {
    f(value);
    f(slope);
  }
  template <class F>
  void for_each_active(F&& f) const {
    f(value);
    f(slope);
  }
};

}

#endif
