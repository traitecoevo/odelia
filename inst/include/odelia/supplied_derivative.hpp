// -*-c++-*-
#ifndef ODELIA_SUPPLIED_DERIVATIVE_HPP_
#define ODELIA_SUPPLIED_DERIVATIVE_HPP_

#include <XAD/XAD.hpp>
#include <XAD/CheckpointCallback.hpp>
#include <odelia/ode_util.hpp>

#include <vector>

namespace odelia {
namespace ode {

// A reverse-sweep callback that supplies a known derivative in place of recorded
// operations. When the sweep reaches it, it reads the adjoint accumulated on the
// output slot and adds dy/dx_i * adjoint to each input slot. Built on XAD's
// CheckpointCallback; construct one through supplied_derivative() below.
template <typename Tape>
class SuppliedDerivative : public xad::CheckpointCallback<Tape> {
  using slot_type = typename Tape::slot_type;
public:
  SuppliedDerivative(slot_type output, std::vector<slot_type> inputs,
               std::vector<double> partials)
    : output_(output), inputs_(std::move(inputs)), partials_(std::move(partials)) {}

  void computeAdjoint(Tape* tape) override {
    const auto ybar = tape->getAndResetOutputAdjoint(output_);
    for (std::size_t i = 0; i < inputs_.size(); ++i) {
      tape->incrementAdjoint(inputs_[i], ybar * partials_[i]);
    }
  }

private:
  slot_type output_;
  std::vector<slot_type> inputs_;
  std::vector<double> partials_;
};

// Make an off-tape result active, carrying its known derivatives into the reverse
// tape. Use this when a value the forward pass needs was computed without recording
// its operations -- a root-find, an optimiser result -- but its derivative is known
// analytically. `y_value` is that value and `partials[i]` is d(y)/d(inputs[i]). `y`
// becomes a fresh tape leaf; a SuppliedDerivative then distributes its adjoint to the
// inputs on the reverse sweep, so downstream code differentiates through `y` without
// the solve being recorded. The tape owns the callback (freed with the tape).
template <typename Tape, typename Active>
Active supplied_derivative(Tape& tape, double y_value,
                     const std::vector<Active*>& inputs,
                     const std::vector<double>& partials) {
  util::check_length(partials.size(), inputs.size());

  Active y = y_value;
  tape.registerInput(y);  // a fresh leaf; downstream operations record its slot

  std::vector<typename Tape::slot_type> input_slots;
  input_slots.reserve(inputs.size());
  for (const auto* x : inputs) {
    input_slots.push_back(x->getSlot());
  }

  auto* edge = new SuppliedDerivative<Tape>(y.getSlot(), std::move(input_slots), partials);
  tape.pushCallback(edge);    // hand ownership to the tape
  tape.insertCallback(edge);  // place the edge at the current tape position
  return y;
}

} // namespace ode
} // namespace odelia

#endif
