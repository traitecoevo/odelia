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
#include <algorithm>
#include <string>
#include <vector>
#include <odelia/ode_callback_system.hpp>

namespace odelia {
namespace rinterface {

namespace detail {

inline void throw_if_domain_sentinel(SEXP out) {
  SEXP msg = Rf_getAttrib(out, Rf_install("odelia_domain_error"));
  if (msg != R_NilValue) {
    util::stop_domain(Rcpp::as<std::string>(msg));
  }
}

// One R callback, called as fn(t, y) or, with parms, as fn(t, y, parms) --
// deSolve's shape, so a function written for it runs unchanged -- with a
// result that may be a list whose first element is the answer (deSolve's
// shape again). The call carries the `odelia_callback` mark that
// domain_error() looks for in sys.calls().
//
// The call is `rhs(t, y, parms)` with `rhs` (or `jac`, `state_valid`) and
// `parms` bound as symbols in a private environment, and t and y placed as
// values. Binding rather than inlining the function means an error raised in
// it reports `rhs(0.14, c(...))` rather than the whole deparsed closure, and
// binding parms means a parms that is itself a symbol or a call reaches the
// function as given rather than being evaluated by the call.
//
// This is the hot path of an R right-hand side, so it does per call only what
// it must: a fresh scalar for t and a fresh vector for y (the callback may keep
// either, so neither can be reused in place), the evaluation, and a copy out.
// The call object and the environment are made once. Results are held with
// PROTECT, not the precious list.
//
// The evaluation is Rcpp_fast_eval: R's unwind-protect with a token allocated
// per call, so that an R error (or domain_error()'s return) unwinds these C++
// frames as an Rcpp::LongjumpException and is resumed at the package boundary.
// A hand-rolled R_UnwindProtect over a cached token and setjmp/longjmp was
// tried in its place and crashed the R process on Windows in three CI runs
// out of three; the per-call token costs nothing measurable.
class Callback {
public:
  Callback(Rcpp::Function fn, Rcpp::RObject parms, const char* name)
      : fn_(fn), parms_(parms), env_(R_NewEnv(R_BaseEnv, FALSE, 0)) {
    SEXP fn_sym = Rf_install(name);
    Rf_defineVar(fn_sym, fn, env_);
    if (parms.isNULL()) {
      call_ = Rf_lang3(fn_sym, R_NilValue, R_NilValue);
    } else {
      SEXP parms_sym = Rf_install("parms");
      Rf_defineVar(parms_sym, parms, env_);
      call_ = Rf_lang4(fn_sym, R_NilValue, R_NilValue, parms_sym);
    }
    Rcpp::Shield<SEXP> mark(Rf_ScalarLogical(1));
    Rf_setAttrib(call_, Rf_install("odelia_callback"), mark);
  }

  // The result, protected by the caller's Shield: a numeric vector or matrix,
  // or whatever the callback returned (checked by the caller).
  SEXP operator()(double t, const std::vector<double>& y) const {
    Rcpp::Shield<SEXP> tt(Rf_ScalarReal(t));
    Rcpp::Shield<SEXP> yy(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(y.size())));
    std::copy(y.begin(), y.end(), REAL(yy));
    SETCADR(call_, tt);
    SETCADDR(call_, yy);
    SEXP out = Rcpp::Rcpp_fast_eval(call_, env_);
    // Leave no dangling reference to this stage's arguments in the call.
    SETCADR(call_, R_NilValue);
    SETCADDR(call_, R_NilValue);
    Rcpp::Shield<SEXP> keep(out);
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
  Rcpp::RObject env_;
  Rcpp::RObject call_;
};

// A REALSXP from a callback result: the result itself, or an integer result
// coerced once. The caller protects what comes back.
inline SEXP as_real(SEXP x) {
  if (TYPEOF(x) == REALSXP) {
    return x;
  }
  if (TYPEOF(x) == INTSXP || TYPEOF(x) == LGLSXP) {
    return Rf_coerceVector(x, REALSXP);
  }
  util::stop("callback must return a numeric vector");
}

} // namespace detail


inline ode::CallbackSystem::rhs_type wrap_rhs(Rcpp::Function rhs, Rcpp::RObject parms) {
  detail::Callback call(rhs, parms, "rhs");
  return [call](double t, const std::vector<double>& y, std::vector<double>& dydt) {
    Rcpp::Shield<SEXP> out(call(t, y));
    Rcpp::Shield<SEXP> rates(detail::as_real(out));
    if (static_cast<size_t>(Rf_xlength(rates)) != y.size()) {
      util::stop("rhs returned " + std::to_string(Rf_xlength(rates)) +
                 " rates, expected " + std::to_string(y.size()));
    }
    std::copy(REAL(rates), REAL(rates) + y.size(), dydt.begin());
  };
}

// jac(t, y) returns an n x n matrix with column j = d f / d y_j; the core wants
// it row-major, J[row * n + col].
inline ode::CallbackSystem::jac_type wrap_jac(Rcpp::Function jac, Rcpp::RObject parms) {
  detail::Callback call(jac, parms, "jac");
  return [call](double t, const std::vector<double>& y,
                const std::vector<double>& /* dydt */, std::vector<double>& J) {
    const size_t n = y.size();
    Rcpp::Shield<SEXP> out(call(t, y));
    Rcpp::Shield<SEXP> m(detail::as_real(out));
    SEXP dim = Rf_getAttrib(m, R_DimSymbol);
    const bool square = TYPEOF(dim) == INTSXP && Rf_xlength(dim) == 2 &&
                        static_cast<size_t>(INTEGER(dim)[0]) == n &&
                        static_cast<size_t>(INTEGER(dim)[1]) == n;
    if (!square) {
      std::string got = TYPEOF(dim) == INTSXP && Rf_xlength(dim) == 2
                            ? std::to_string(INTEGER(dim)[0]) + " x " + std::to_string(INTEGER(dim)[1]) + " matrix"
                            : "vector of length " + std::to_string(Rf_xlength(m));
      util::stop("jac returned a " + got + ", expected " + std::to_string(n) +
                 " x " + std::to_string(n));
    }
    const double* v = REAL(m); // column-major
    for (size_t row = 0; row < n; ++row) {
      for (size_t col = 0; col < n; ++col) {
        J[row * n + col] = v[row + col * n];
      }
    }
  };
}

inline ode::CallbackSystem::valid_type wrap_valid(Rcpp::Function valid, Rcpp::RObject parms) {
  detail::Callback call(valid, parms, "state_valid");
  return [call](double t, const std::vector<double>& y) {
    Rcpp::Shield<SEXP> out(call(t, y));
    if (TYPEOF(out) != LGLSXP || Rf_xlength(out) != 1 || LOGICAL(out)[0] == NA_LOGICAL) {
      util::stop("state_valid must return a single TRUE or FALSE");
    }
    return LOGICAL(out)[0] == 1;
  };
}

inline ode::CallbackSystem make_r_system(Rcpp::Function rhs,
                                         Rcpp::Nullable<Rcpp::Function> jac,
                                         Rcpp::Nullable<Rcpp::Function> valid,
                                         Rcpp::RObject parms,
                                         std::vector<double> y0, double t0,
                                         bool autonomous, double jac_fd_step,
                                         double jac_fd_floor) {
  ode::CallbackSystem::jac_type j;
  if (jac.isNotNull()) {
    j = wrap_jac(Rcpp::Function(jac.get()), parms);
  }
  ode::CallbackSystem::valid_type v;
  if (valid.isNotNull()) {
    v = wrap_valid(Rcpp::Function(valid.get()), parms);
  }
  return ode::CallbackSystem(wrap_rhs(rhs, parms), std::move(y0), t0, j, v,
                             autonomous, jac_fd_step, jac_fd_floor);
}

} // namespace rinterface
} // namespace odelia

#endif
