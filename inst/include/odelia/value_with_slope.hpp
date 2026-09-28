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
// it. plant's competition path carries d(value)/d(height); phylloptim's operating
// point carries d(value)/d(collar potential). Naming the variable here would make
// one of them wrong.
//
// ⚠️ IT LIVES HERE BECAUSE OF for_each_active, NOT BECAUSE IT IS SHARED. visit_active
// passes over any shape it does not open, without refusing it, and it does not open
// an aggregate of two scalars -- so a pair that does not say what it holds loses
// both members from the walk, silently, and active_system::release is what reports
// the miss. That obligation is this library's, so a model defining the pair itself
// is a model carrying this library's problem.
//
// No includes of its own, so a header on any path can take it without taking
// anything else with it.
template <typename T>
struct with_slope {
  T value;
  T slope;

  // Both, because visit_active dispatches on whether the call compiles: handed a
  // const object it drops an arm that is non-const and passes over the whole
  // shape in silence. The rewinding forms in implicit_node.hpp take
  // `const Inputs&...` -- they read and clear through the tape by slot, so const
  // is what they want -- and without the const overload here a pair handed to one
  // of them contributes NO rows, which arrives as an exact zero in a column.
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
