#ifndef EMHAWKES_R_PRECISION_H
#define EMHAWKES_R_PRECISION_H

#include <Rcpp.h>
#include <limits>

namespace emhawkes {
// The compiler may offer long double even when R was configured without it.
// Query R once per DLL, before entering an event loop or an RNG scope.
inline bool use_r_long_double() {
    static const bool enabled = []() {
        Rcpp::Function capabilities = Rcpp::Environment::base_env()["capabilities"];
        return Rcpp::as<bool>(capabilities("long.double"));
    }();
    return enabled;
}

// R materializes vector products before summing. With double accumulation,
// prevent contraction of a product and its following addition into an FMA.
template <typename Sum>
inline double rounded_value(double value) {
    if (std::numeric_limits<Sum>::digits == std::numeric_limits<double>::digits) {
        volatile double rounded = value;
        return rounded;
    }
    return value;
}

template <typename Sum>
inline double rounded_product(double x, double y) {
    return rounded_value<Sum>(x * y);
}
}
#endif
