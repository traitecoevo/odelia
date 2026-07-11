/* An example of odelia::ode::supplied_derivative, compiled on demand by
 * test-ad-supplied-derivative.R (sourceCpp).
 *
 * A value computed OFF the tape (a Newton root-find) is made active by supplying its
 * implicit-function-theorem partials, then a downstream expression is differentiated
 * through it without the solve being recorded. The test compares against the
 * closed-form derivative and finite differences.
 */

// [[Rcpp::depends(Rcpp, odelia)]]
// [[Rcpp::plugins(cpp20)]]

#include <Rcpp.h>
#include <XAD/XAD.hpp>
#include <odelia/supplied_derivative.hpp>

#include <cmath>

using namespace Rcpp;
using namespace odelia;

// Off-tape Newton solve of x = cos(x) + a; returns the root x*(a). Deliberately
// a transcendental with no closed form for x, so differentiating it *through* the
// iterations is exactly what the supplied derivative lets us avoid.
static double solve_root(double a) {
  double x = a;
  for (int i = 0; i < 100; ++i) {
    const double f = x - std::cos(x) - a;
    const double dx = f / (1.0 + std::sin(x)); // f'(x) = 1 + sin(x)
    x -= dx;
    if (std::abs(dx) < 1e-15) break;
  }
  return x;
}

// [[Rcpp::export]]
Rcpp::List supplied_derivative_demo(double a) {
  using ad = xad::adj<double>;
  using ad_type = ad::active_type;

  ad::tape_type tape;
  ad_type A = a;
  tape.registerInput(A);
  tape.newRecording();

  // Root-find on the value only; its IFT sensitivity dx/da = 1 / (1 + sin(x*))
  // follows from F(x, a) = x - cos(x) - a = 0, so dx/da = -F_a / F_x.
  const double xv = solve_root(xad::value(A));
  const double dxda = 1.0 / (1.0 + std::sin(xv));

  // Inject x as an active leaf carrying dx/dA, then differentiate g = x^2.
  ad_type x = ode::supplied_derivative(tape, xv, std::vector<ad_type*>{&A},
                                 std::vector<double>{dxda});
  ad_type g = x * x;

  tape.registerOutput(g);
  xad::derivative(g) = 1.0;
  tape.computeAdjoints();

  return Rcpp::List::create(
    Rcpp::Named("x") = xv,
    Rcpp::Named("g") = xad::value(g),
    Rcpp::Named("dg_da") = xad::derivative(A),
    Rcpp::Named("dg_da_analytic") = 2.0 * xv * dxda);
}
