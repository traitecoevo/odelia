#ifndef ODELIA_RCPP_INTERFACE_HELPERS_HPP_
#define ODELIA_RCPP_INTERFACE_HELPERS_HPP_

#include <Rcpp.h>
#include <odelia/drivers.hpp>

namespace odelia {

inline Rcpp::XPtr<drivers::Drivers> get_Drivers(SEXP xp) {
  return Rcpp::XPtr<drivers::Drivers>(xp);
}

}  // namespace odelia

#endif
