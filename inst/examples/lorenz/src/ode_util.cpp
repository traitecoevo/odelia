#include <odelia/ode_util.hpp>
#include <Rcpp.h>

namespace odelia
{
  namespace util
  {

    void check_length(size_t received, size_t expected)
    {
      if (expected != received)
      {
        Rcpp::stop("Incorrect length input; expected " +
                   std::to_string(expected) + ", received " +
                   std::to_string(received));
      }
    }

    bool is_finite(double x)
    {
      return R_FINITE(x);
    }

    void stop(const std::string &msg)
    {
      Rcpp::stop(msg);
    }

    void warning(const std::string &msg)
    {
      Rcpp::warning(msg);
    }

  }
}

namespace Rcpp {
template <> SEXP wrap(const odelia::util::index& x) {
  return Rcpp::wrap(odelia::util::base_0_to_1<size_t, int>(x.x));
}
template <> odelia::util::index as(SEXP x) {
  const int ix(Rcpp::as<int>(x));
  if (ix <= 0) {
    Rcpp::stop("Invalid value for index (must be >= 1)");
  }
  return odelia::util::base_1_to_0<int, size_t>(ix);
}
template <> SEXP wrap(const std::vector<odelia::util::index>& x) {
  Rcpp::IntegerVector ret(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[static_cast<int>(i)] = odelia::util::base_0_to_1<size_t, int>(x[i].x);
  }
  return Rcpp::wrap(ret);
}
}
