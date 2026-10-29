// -*-c++-*-
#ifndef ODELIA_ADJOINT_SCALAR_HPP_
#define ODELIA_ADJOINT_SCALAR_HPP_

#include <XAD/XAD.hpp>
#include <concepts>
#include <cstddef>
#include <type_traits>

// The reverse-mode scalar and the walk over the active values a class holds,
// apart from the System contract in ode_interface.hpp because a consumer that
// records partials on a tape needs these and no solver.

namespace odelia {
namespace ode {

// The scalar a reverse-mode (adjoint) pass runs on: one adjoint layer above T.
// This is what a System is rebound to for a gradient, and the type its output
// adjoints and recorded values have. T is a parameter so the layer can sit on
// another active scalar; at the default it sits on double, which is the scalar
// an ordinary solve is differentiated from.
//
// One derivative per slot, deliberately. `xad::adj` takes a second, defaulted
// template argument for the derivative width, so `xad::adj<T, 3>` would sweep
// three outputs in one walk. It does not pay: widening multiplies the bytes the
// derivative array occupies, and the cache gives back more than the shared
// traversal saves. Every width agrees bit for bit, so what is ruled out is the
// cost and never the answer.
template <typename T = double>
using adjoint_scalar = typename xad::adj<T>::active_type;

// The tape a reverse pass records on, reached from the scalar rather than named
// beside it. Naming a tape independently of its scalar is right at one
// derivative width and silently wrong at another -- both spellings compile, and
// the one that queries the wrong tape reports no tape active and stops nothing.
template <typename T = double>
using adjoint_tape = typename adjoint_scalar<T>::tape_type;

// Every active value among the members handed in, whatever shape they arrive in:
// a scalar, a container of them, a pointer to one, or an object that answers
// for_each_active. This skips a member the visitor cannot be called with, so a
// class lists what it holds rather than deciding which of them carry a
// derivative -- which is the judgment that makes forgetting one likely.
//
// The skip is what carries a passive member: the visitor takes the active scalar,
// so a double or a bool matches no arm at all and this leaves it alone. A class
// holding doubles only declares no for_each_active, and is walked past the same
// way.
//
// ⚠️ A MEMBER IN A SHAPE NOT LISTED HERE IS SKIPPED, NOT REFUSED. An active scalar
// inside something this does not open is passed over in silence, and the count a
// caller checks afterwards is the only thing that says so. This cannot refuse the
// shape instead: from here a class holding active values and a class holding none
// look the same, and callers hand over both. What reports the miss is
// adjoint_system::release, which counts the slots the walk did not reach.
//
// This does not open an aggregate of two scalars, so a type standing in for a
// std::pair on a recorded path declares for_each_active, as hermite_spline and
// with_slope do.
template <class F, class T>
void visit_active(F& f, T& x) {
  // A type that declares for_each_active but cannot be walked const would stop
  // matching the first arm and fall through in silence, and the rewinding forms
  // in implicit_node.hpp hand their inputs over const.
  static_assert(
      !(std::is_const_v<T> &&
        requires(std::remove_const_t<T>& mutable_x) { mutable_x.for_each_active(f); } &&
        !requires { x.for_each_active(f); }),
      "visit_active: this type declares for_each_active but cannot be walked "
      "through a const reference, so every active scalar it holds would be "
      "skipped in silence. Add a const overload of for_each_active.");
  if constexpr (requires { x.for_each_active(f); }) {
    x.for_each_active(f);
  } else if constexpr (requires { x.begin(); x.end(); }) {
    for (auto& element : x) { visit_active(f, element); }
  } else if constexpr (requires { *x; x == nullptr; }) {
    if (x != nullptr) { visit_active(f, *x); }
  } else if constexpr (std::invocable<F&, T&>) {
    f(x);
  }
}

template <class F, class A, class B, class... Rest>
void visit_active(F& f, A& a, B& b, Rest&... rest) {
  visit_active(f, a);
  visit_active(f, b, rest...);
}

// How many of these hold a tape slot -- which is what the walk above can reach
// and no more. A test that knows its own input list asserts this beside the call
// it is testing; production code has nothing to compare it against, because the
// only other count of the same thing is this walk.
//
// It says nothing about a value the list LEAVES OUT, or a member a type's
// for_each_active omits: neither is walked, so neither is counted.
template <class S, class... Inputs>
std::size_t count_active_slots(const Inputs&... inputs) {
  std::size_t n = 0;
  auto count = [&](const S& x) {
    if (x.getSlot() != S::tape_type::INVALID_SLOT) {
      ++n;
    }
  };
  visit_active(count, inputs...);
  return n;
}

}
}

#endif
