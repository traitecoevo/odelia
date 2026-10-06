#ifndef ODELIA_RCPP_INTERFACE_HELPERS_HPP_
#define ODELIA_RCPP_INTERFACE_HELPERS_HPP_

#include <Rcpp.h>
#include <string>
#include <odelia/drivers.hpp>
#include <odelia/ode_solver_internal.hpp>

namespace odelia {

inline Rcpp::XPtr<drivers::Drivers> get_Drivers(SEXP xp) {
  return Rcpp::XPtr<drivers::Drivers>(xp);
}

// Map an R-facing method string to the solver Method enum.
inline ode::Method parse_method(const std::string& method) {
  if (method == "rodas" || method == "implicit") {
    return ode::Method::rodas;
  }
  if (method == "rkck" || method == "rk45" || method == "explicit") {
    return ode::Method::rkck;
  }
  if (method == "dopri" || method == "dopri5" || method == "ode45" || method == "rk45dp7") {
    return ode::Method::dopri;
  }
  Rcpp::stop("Unknown method '" + method + "'. Use 'dopri', 'rkck' or 'rodas'.");
}

}  // namespace odelia

#endif
