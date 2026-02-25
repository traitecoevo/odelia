#include <Rcpp.h>
#include <examples/parabola_system.hpp>

// [[Rcpp::export]]
Rcpp::List parabola_eval(Rcpp::NumericVector params) {
    std::vector<double> p(params.begin(), params.end());
    auto result = parabola_with_gradient(p);
    return Rcpp::List::create(
        Rcpp::Named("value") = result.first,
        Rcpp::Named("gradient") = Rcpp::wrap(result.second)
    );
}
