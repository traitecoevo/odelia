// -*-c++-*-
#ifndef ODELIA_TANGENT_HPP_
#define ODELIA_TANGENT_HPP_

#include <XAD/XAD.hpp>

// The forward-mode vocabulary, apart from the reverse-mode one because a
// consumer can want this and not that: a package with no tape reaches for the
// tangent scalar and must not be handed a name for a tape it does not have.

namespace odelia {
namespace ode {

// The scalar a forward-mode (tangent) pass runs on: one tangent layer above T,
// carrying a directional derivative and no tape. Nest it for a curvature -- a
// tangent of a tangent -- which is the only second-order scalar this family
// uses.
//
// A tangent above an adjoint is refused. `BinaryExpr` stores its operands and
// its cached value BY VALUE, and `value()` and `derivative()` return the scalar
// by value, so at an active inner scalar every one of those copies is a recorded
// statement, and the count grows faster than the expression's depth. A slope
// wanted AT an adjoint scalar is taken through the kernel at double and handed
// over as supplied rows -- `record_with_derivatives`, one statement whatever the
// row count. At any nested scalar, flatten an expression into named
// intermediates rather than fusing it into one.
template <typename T>
struct tangent_over {
  static_assert(!xad::ExprTraits<T>::isReverse,
                "a tangent above an adjoint records every operand copy; take "
                "the slope through the kernel at double and supply the row");
  using type = typename xad::fwd<T>::active_type;
};
template <typename T = double>
using tangent_scalar = typename tangent_over<T>::type;

// A scalar that carries a direction rather than an adjoint accumulator.
template <typename S>
concept CarriesDirection = xad::ExprTraits<S>::isForward;

// The direction a forward pass carries in, and the derivative it carries out.
// Everything about them is a pass-through except which scalars they accept.
//
// That is the whole reason they are named. One accessor spells this AND an
// adjoint scalar's accumulator, so on an adjoint scalar the same statement seeds
// a slot no forward pass reads, or reads one no sweep has written -- and neither
// raises anything, because both are the accessor doing exactly what it says.
template <typename S>
  requires CarriesDirection<S>
void seed_direction(S& x, double direction) {
  xad::derivative(x) = direction;
}

// Returned by value, and declared rather than deduced. By value because the
// argument is routinely a temporary -- the kernel whose derivative is being read
// -- and a reference into one is a reference a caller can outlive. It costs
// nothing to copy: a tangent carries no tape, so a copy of one is bytes, where a
// copy of an adjoint value would record an operation.
template <typename S>
  requires CarriesDirection<S>
typename S::derivative_type derivative_along(const S& x) {
  return xad::derivative(x);
}

}
}

#endif
