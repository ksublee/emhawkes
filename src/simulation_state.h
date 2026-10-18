#ifndef EMHAWKES_SIMULATION_STATE_H
#define EMHAWKES_SIMULATION_STATE_H

#include <Rcpp.h>
#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <string>
#include <vector>

namespace {
// Simulation may be nested inside an R warning handler or callback. Always
// synchronize explicitly rather than relying on Rcpp's global scope counter.
class SimulationRNGScope {
public:
    SimulationRNGScope() { GetRNGstate(); }
    ~SimulationRNGScope() { PutRNGstate(); }
    SimulationRNGScope(const SimulationRNGScope&) = delete;
    SimulationRNGScope& operator=(const SimulationRNGScope&) = delete;
};

// Direct Rcpp::warning() can long-jump over destructors under warn=2. Protected
// R evaluation preserves conditions and unwinds C++ state. Handlers may draw
// random numbers or simulate another path, so publish and reload RNG state.
inline void simulation_warning(const std::string& message) {
    PutRNGstate();
    try {
        Rcpp::Function warn = Rcpp::Environment::base_env()["warning"];
        warn(message, Rcpp::Named("call.") = false);
        GetRNGstate();
    } catch (...) {
        GetRNGstate();
        throw;
    }
}

// R materializes double vector products before sum(). On platforms where
// long double is double, prevent contraction into a fused multiply-add.
#if LDBL_MANT_DIG == DBL_MANT_DIG
inline double rounded_product(double x, double y) {
    volatile double value = x * y;
    return value;
}
#else
// Keep the original double multiplication on extended-precision platforms.
#define rounded_product(x, y) ((x) * (y))
#endif

// R's rowSums accumulates in long double, then rounds once to double.
double row_sum(const std::vector<double>& state, int row, int d, int k) {
    long double sum = 0.0L;
    for (int j = 0; j < k; ++j) sum += state[row + d * j];
    return static_cast<double>(sum);
}

typedef std::vector<std::vector<std::vector<int> > > RateGroups;

// Preserve the usual expression exactly; only subnormal intermediate values
// use the continuous small-decay limit instead of quantizing waiting times.
double excitation_wait(double beta, double log_u, double total) {
    const double product = beta * log_u;
    const double scaled = product / total;
    if (total <= 0.0) return R_PosInf;
    if (std::abs(product) < DBL_MIN || std::abs(scaled) < DBL_MIN) {
        const double first_order = -log_u / total;
        const double q = beta * first_order;
        if (!R_FINITE(q) || q >= 1.0) return R_PosInf;
        return first_order * (q > 0.0 ? -std::log1p(-q) / q : 1.0);
    }
    return scaled > -1.0 ? -std::log1p(scaled) / beta : R_PosInf;
}

bool valid_state(const std::vector<double>& state, int d, int k,
    const Rcpp::NumericMatrix& beta, const RateGroups& groups) {
    for (int i = 0; i < d; ++i) {
        bool signed_row = false;
        long double sum = 0.0L, scale = 0.0L;
        for (int j = 0; j < k; ++j) {
            const double value = state[i + d * j];
            if (!R_FINITE(value)) return false;
            signed_row = signed_row || value < 0.0;
            sum += value;
            scale += std::abs(value);
        }
        if (!signed_row) continue;
        if (groups[i].size() == 1) {
            // Match the separate rounded sum()/sum(abs()) in R validation.
            if (!R_FINITE(static_cast<double>(scale)) ||
                static_cast<double>(sum) < -32.0 * DBL_EPSILON *
                    std::max(1.0, static_cast<double>(scale))) return false;
        } else {
            scale = 0.0L;
            for (int j = 0; j < k; ++j)
                scale += std::abs(rounded_product(state[i + d * j], beta(i, j)));
            if (!R_FINITE(static_cast<double>(scale))) return false;
            const double tolerance = -32.0 * DBL_EPSILON *
                std::max(1.0, static_cast<double>(scale));
            long double cumulative = 0.0L;
            for (const auto& group : groups[i]) {
                long double weight = 0.0L;
                for (int j : group) weight += rounded_product(state[i + d * j], beta(i, j));
                // R rounds each group sum before cumsum's extended accumulation.
                cumulative += static_cast<double>(weight);
                if (!R_FINITE(static_cast<double>(cumulative)) ||
                    static_cast<double>(cumulative) < tolerance) return false;
            }
        }
    }
    return true;
}
}


#endif
