#include <Rcpp.h>
#include <climits>
#include <cmath>
#include "jacobian_call.h"

// Retain event-score projections without interleaving baseline and impact calls.
// Small dimensions keep the dense path; the initial Jacobian remains separate.
// [[Rcpp::export(rng = false)]]
Rcpp::RObject hawkes_jacobian_path_cpp(SEXP callback, Rcpp::List arguments,
    int dimensions, Rcpp::NumericVector dt, Rcpp::IntegerVector type,
    Rcpp::CharacterVector parameter_names) {
    const int parameters = parameter_names.size();
    if (dimensions < 1 || dt.size() < 1 || dt.size() > INT_MAX ||
        type.size() != dt.size() || (parameters && dimensions > INT_MAX / parameters))
        Rcpp::stop("Invalid native baseline Jacobian dimensions.");
    const int size = static_cast<int>(dt.size());
    for (int n = 0; n < size; ++n) {
        if ((n & 65535) == 0) Rcpp::checkUserInterrupt();
        if (!R_FINITE(dt[n]) || dt[n] < 0 || (n && (type[n] < 1 || type[n] > dimensions)))
            Rcpp::stop("Invalid native baseline Jacobian observations.");
    }
    emhawkes::JacobianCall call(callback, arguments, dimensions, parameters,
                                parameter_names, true, false, "jacobian$mu");
    const bool compact = dimensions > 2 &&
        static_cast<long long>(size) * (2LL * dimensions - 5) > 2LL * dimensions;
    Rcpp::NumericMatrix path(compact ? parameters : dimensions * parameters, size), dummy(0, 0);
    Rcpp::NumericMatrix initial(compact ? dimensions : 0, compact ? parameters : 0);
    Rcpp::NumericMatrix sum(compact ? parameters : 0, compact ? size : 0);
    Rcpp::IntegerMatrix exponent(compact ? parameters : 0, compact ? size : 0);
    std::vector<double> state;
    for (int n = 0; n < size; ++n) {
        if ((n & 1023) == 0) Rcpp::checkUserInterrupt();
        call.evaluate(n + 1, dummy, state);
        for (int p = 0; p < parameters; ++p) {
            if (!compact) {
                for (int i = 0; i < dimensions; ++i) path(i + dimensions * p, n) = call.at(i, p);
            } else if (n == 0) {
                for (int i = 0; i < dimensions; ++i) initial(i, p) = call.at(i, p);
            } else {
                path(p, n) = call.at(type[n] - 1, p);
                long double total = 0.0L;
                for (int i = 0; i < dimensions; ++i) total += call.at(i, p);
                int power = 0;
                // On platforms where long double has double range, retry an
                // overflowing sum at a common binary scale. Ordinary values
                // need neither a maximum scan nor per-dimension rescaling.
                if (!std::isfinite(total)) {
                    double largest = 0.0;
                    for (int i = 0; i < dimensions; ++i)
                        largest = std::max(largest, std::abs(call.at(i, p)));
                    std::frexp(largest, &power);
                    total = 0.0L;
                    for (int i = 0; i < dimensions; ++i)
                        total += std::ldexp(static_cast<long double>(call.at(i, p)), -power);
                }
                // A scaled sum avoids overflow and loss of subnormal inputs.
                // The event loop restores its exponent after multiplying by dt.
                int sum_power = 0;
                sum(p, n) = static_cast<double>(std::frexp(total, &sum_power));
                exponent(p, n) = power + sum_power;
            }
        }
    }
    if (compact) return Rcpp::List::create(Rcpp::Named("initial") = initial,
        Rcpp::Named("event") = path, Rcpp::Named("sum") = sum, Rcpp::Named("exponent") = exponent);
    return path;
}
