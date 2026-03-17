#ifndef ODELIA_ODE_FIT_HPP_
#define ODELIA_ODE_FIT_HPP_

#include <XAD/XAD.hpp>
#include <odelia/ode_solver.hpp>
#include <optional>
#include <vector>
#include <cmath>

namespace odelia {
namespace ode {

// Loss function
template<typename T>
T sum_of_squares(const std::vector<std::vector<T>>& obs, 
                 const std::vector<std::vector<double>>& targets) {
    T loss(0.0);
    for (size_t i = 0; i < obs.size(); ++i) {
        for (size_t j = 0; j < obs[i].size(); ++j) {
            T diff = obs[i][j] - targets[i][j]; 
            loss += diff * diff;
        }
    }
    return loss;
}

// Compute gradient w.r.t. initial conditions and/or parameters
template<typename Solver>
std::pair<double, std::vector<double>> compute_gradient(
    Solver& solver,
    std::optional<std::vector<double>> ic = std::nullopt,
    std::optional<std::vector<double>> params = std::nullopt
) {
    using ad = xad::adj<double>;
    using ad_type = ad::active_type;
    
    // At least one must be provided
    if (!ic && !params) {
        util::stop("Must provide at least one of 'ic' or 'params'");
    }
    
    // Create tape
    ad::tape_type tape;

    // Collect input pointers for gradient extraction
    std::vector<ad_type*> inputs;

        constexpr double nonfinite_penalty = 1e100;
        auto penalty_result = [&inputs, nonfinite_penalty]() {
            return std::make_pair(nonfinite_penalty, std::vector<double>(inputs.size(), 0.0));
        };

    auto& system = solver.get_system_ref();

    // Set params and register BEFORE newRecording (XAD requirement)
    if (params) {
        auto refs = system.set_params(tape, params->begin());
        inputs.insert(inputs.end(), refs.begin(), refs.end());
    }

    // Set initial conditions and register BEFORE newRecording (XAD requirement)
    if (ic) {
        auto refs = system.set_initial_state(tape, ic->begin(), solver.fit_times()[0]);
        inputs.insert(inputs.end(), refs.begin(), refs.end());
    }

    // Start recording AFTER registering inputs (XAD requirement)
    tape.newRecording();

    // Reset AFTER recording starts so the copy operation is recorded as a dependency
    // This ensures compute_rates() in reset() uses the registered AD types
    solver.reset();
    
    // Forward pass - returns states at observation times
    auto obs = solver.advance_target();

        // If trajectory values become non-finite, return a large finite loss so
        // optimizers can continue searching without aborting on NaN/Inf.
        for (const auto& row : obs) {
            for (const auto& value : row) {
                if (!std::isfinite(xad::value(value))) {
                    return penalty_result();
                }
            }
        }
    
    // Compute loss
    ad_type loss = sum_of_squares(obs, solver.targets());

        if (!std::isfinite(xad::value(loss))) {
            return penalty_result();
        }
    
    // Backpropagate
    tape.registerOutput(loss);
    xad::derivative(loss) = 1.0;
    tape.computeAdjoints();
    
    // Extract gradients from registered inputs
    std::vector<double> gradient(inputs.size());
    for (size_t i = 0; i < inputs.size(); ++i) {
                const double g = xad::derivative(*inputs[i]);
                if (!std::isfinite(g)) {
                    return penalty_result();
                }
                gradient[i] = g;
    }
    
    return {xad::value(loss), gradient};
}

} // namespace ode
} // namespace odelia

#endif
