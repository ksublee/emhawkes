#include <Rcpp.h>
#include <climits>
#include "jacobian_call.h"

// Pack observed-baseline total Jacobians in column-major order at each index.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hawkes_jacobian_path_cpp(SEXP callback, Rcpp::List arguments,
    int dimensions, int size, Rcpp::CharacterVector parameter_names) {
    const int parameters = parameter_names.size();
    if (dimensions < 1 || size < 1 || (parameters && dimensions > INT_MAX / parameters))
        Rcpp::stop("Invalid native baseline Jacobian dimensions.");
    emhawkes::JacobianCall call(callback, arguments, dimensions, parameters,
                                parameter_names, true, false, "jacobian$mu");
    Rcpp::NumericMatrix path(dimensions * parameters, size), dummy(0, 0);
    std::vector<double> state;
    for (int n = 0; n < size; ++n) {
        if ((n & 1023) == 0) Rcpp::checkUserInterrupt();
        call.evaluate(n + 1, dummy, state);
        for (int p = 0; p < parameters; ++p)
            for (int i = 0; i < dimensions; ++i) path(i + dimensions * p, n) = call.at(i, p);
    }
    return path;
}
