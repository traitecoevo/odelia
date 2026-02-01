#ifndef PARABOLA_SYSTEM_HPP
#define PARABOLA_SYSTEM_HPP

#include <XAD/XAD.hpp>
#include <vector>
#include <utility>

// Pure objective - templated for double or AReal<double>
template <typename T>
T parabola_objective(T x, T y) {
    return -((x - 2.0) * (x - 2.0)) - ((y - 3.0) * (y - 3.0));
}

// AD wrapper - encapsulates tape logic
inline std::pair<double, std::vector<double>>
parabola_with_gradient(const std::vector<double>& params) {
    using ad = xad::adj<double>;
    ad::tape_type tape;
    
    ad::active_type x = params[0], y = params[1];
    tape.registerInput(x);
    tape.registerInput(y);
    tape.newRecording();
    
    ad::active_type f = parabola_objective(x, y);
    
    tape.registerOutput(f);
    xad::derivative(f) = 1.0;
    tape.computeAdjoints();
    
    return {xad::value(f), {xad::derivative(x), xad::derivative(y)}};
}

#endif
