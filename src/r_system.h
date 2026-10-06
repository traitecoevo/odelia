// -*-c++-*-
#ifndef ODELIA_R_SYSTEM_HPP_
#define ODELIA_R_SYSTEM_HPP_

// The R adapter for ode::CallbackSystem (#62): R closures wrapped as the
// std::functions the core system takes. This is the only place R and the
// callback system meet, and it lives in src/ -- the header core stays R-free.
//
// Errors. An R error inside a callback is a bug, and must reach R unchanged:
// class, message and call intact. Rcpp evaluates the callback under R's
// unwind-protect, so the error unwinds these C++ frames as a token and is
// resumed at the package boundary; nothing here needs to see it.
//
// A *domain* refusal -- "this state is impossible, take a smaller step" -- is
// not an error and does not travel as one, because every mechanism that would
// let C++ classify an R condition costs more than the callback itself
// (R-level tryCatch around each call: 2.6x on the Lorenz benchmark;
// R_tryCatchError, which is R-level tryCatch underneath: 7x). Instead the call
// into the callback is marked -- the LANGSXP carries an `odelia_callback`
// attribute, which R's sys.calls() preserves -- and domain_error() finds that
// frame and *returns* from it with a sentinel value, an NA carrying the message
// as an attribute, which the adapter below throws as util::DomainError so the
// adaptive loop rejects the step (#55). The return is a genuine non-local exit
// from any depth below the callback, done by evaluating `return(sentinel)` in
// the frame from C (odelia_return_from): through Rf_eval directly, not R's
// eval(), which would catch the return itself. The fast path pays one
// attribute on the call and nothing else.

#include <Rcpp.h>
#include <string>
#include <vector>
#include <odelia/ode_callback_system.hpp>

namespace odelia {
namespace rinterface {

namespace detail {

inline void throw_if_domain_sentinel(const Rcpp::RObject& out) {
  if (out.hasAttribute("odelia_domain_error")) {
    util::stop_domain(Rcpp::as<std::string>(out.attr("odelia_domain_error")));
  }
}

// One R callback, called as fn(t, y) or, with parms, as fn(t, y, parms) --
// deSolve's shape, so a function written for it runs unchanged -- with a
// result that may be a list whose first element is the answer (deSolve's
// shape again). The call object is built once and its arguments replaced on
// each evaluation: an allocation per stage would be most of what an R
// right-hand side costs. The call carries the `odelia_callback` mark that
// domain_error() looks for in sys.calls().
class Callback {
public:
  Callback(Rcpp::Function fn, Rcpp::RObject parms)
      : fn_(fn), parms_(parms),
        call_(parms.isNULL() ? Rf_lang3(fn, R_NilValue, R_NilValue)
                             : Rf_lang4(fn, R_NilValue, R_NilValue, parms)) {
    Rcpp::Shield<SEXP> mark(Rf_ScalarLogical(1));
    Rf_setAttrib(call_, Rf_install("odelia_callback"), mark);
  }

  Rcpp::RObject operator()(double t, const std::vector<double>& y) const {
    Rcpp::Shield<SEXP> tt(Rf_ScalarReal(t));
    Rcpp::Shield<SEXP> yy(Rcpp::wrap(y));
    SETCADR(call_, tt);
    SETCADDR(call_, yy);
    Rcpp::RObject out(Rcpp::Rcpp_fast_eval(call_, R_GlobalEnv));
    // Leave no dangling reference to this stage's arguments in the call.
    SETCADR(call_, R_NilValue);
    SETCADDR(call_, R_NilValue);
    throw_if_domain_sentinel(out);
    if (TYPEOF(out) == VECSXP) {
      if (Rf_xlength(out) < 1) {
        util::stop("callback returned an empty list");
      }
      out = VECTOR_ELT(out, 0);
      throw_if_domain_sentinel(out);
    }
    return out;
  }

private:
  Rcpp::Function fn_;
  Rcpp::RObject parms_;
  Rcpp::RObject call_;
};

} // namespace detail

inline ode::CallbackSystem::rhs_type wrap_rhs(Rcpp::Function rhs, Rcpp::RObject parms) {
  detail::Callback call(rhs, parms);
  return [call](double t, const std::vector<double>& y, std::vector<double>& dydt) {
    Rcpp::NumericVector rates(call(t, y));
    if (static_cast<size_t>(rates.size()) != y.size()) {
      util::stop("rhs returned " + std::to_string(rates.size()) +
                 " rates, expected " + std::to_string(y.size()));
    }
    for (size_t i = 0; i < y.size(); ++i) {
      dydt[i] = rates[i];
    }
  };
}

// jac(t, y) returns an n x n matrix with column j = d f / d y_j; the core wants
// it row-major, J[row * n + col].
inline ode::CallbackSystem::jac_type wrap_jac(Rcpp::Function jac, Rcpp::RObject parms) {
  detail::Callback call(jac, parms);
  return [call](double t, const std::vector<double>& y,
                const std::vector<double>& /* dydt */, std::vector<double>& J) {
    const size_t n = y.size();
    Rcpp::NumericMatrix m(call(t, y));
    if (static_cast<size_t>(m.nrow()) != n || static_cast<size_t>(m.ncol()) != n) {
      util::stop("jac returned a " + std::to_string(m.nrow()) + " x " +
                 std::to_string(m.ncol()) + " matrix, expected " +
                 std::to_string(n) + " x " + std::to_string(n));
    }
    for (size_t row = 0; row < n; ++row) {
      for (size_t col = 0; col < n; ++col) {
        J[row * n + col] = m(row, col);
      }
    }
  };
}

inline ode::CallbackSystem::valid_type wrap_valid(Rcpp::Function valid, Rcpp::RObject parms) {
  detail::Callback call(valid, parms);
  return [call](double t, const std::vector<double>& y) {
    Rcpp::LogicalVector ok(call(t, y));
    if (ok.size() != 1 || Rcpp::LogicalVector::is_na(ok[0])) {
      util::stop("state_valid must return a single TRUE or FALSE");
    }
    return ok[0] == TRUE;
  };
}

inline ode::CallbackSystem make_r_system(Rcpp::Function rhs,
                                         Rcpp::Nullable<Rcpp::Function> jac,
                                         Rcpp::Nullable<Rcpp::Function> valid,
                                         Rcpp::RObject parms,
                                         std::vector<double> y0, double t0,
                                         bool autonomous, double jac_fd_step) {
  ode::CallbackSystem::jac_type j;
  if (jac.isNotNull()) {
    j = wrap_jac(Rcpp::Function(jac.get()), parms);
  }
  ode::CallbackSystem::valid_type v;
  if (valid.isNotNull()) {
    v = wrap_valid(Rcpp::Function(valid.get()), parms);
  }
  return ode::CallbackSystem(wrap_rhs(rhs, parms), std::move(y0), t0, j, v,
                             autonomous, jac_fd_step);
}

} // namespace rinterface
} // namespace odelia

#endif
