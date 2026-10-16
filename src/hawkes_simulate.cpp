#include <Rcpp.h>
#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <vector>

namespace {
// R's rowSums accumulates in long double, then rounds once to double.
double row_sum(const std::vector<double>& state, int row, int d, int k) {
    long double sum = 0.0L;
    for (int j = 0; j < k; ++j) sum += state[row + d * j];
    return static_cast<double>(sum);
}

bool valid_state(const std::vector<double>& state, int d, int k) {
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
        // Match the separate rounded sum()/sum(abs()) used in R validation.
        if (signed_row && (!R_FINITE(static_cast<double>(scale)) ||
            static_cast<double>(sum) < -32.0 * DBL_EPSILON *
                std::max(1.0, static_cast<double>(scale)))) return false;
    }
    return true;
}
}

// Fresh output/count buffers are prepared privately by simulate_hawkes().
// No user input matrix is mutated. Support is decided in R before this call.
// The explicit scope below avoids even creating .Random.seed for size=1.
// [[Rcpp::export(rng = false)]]
Rcpp::List hawkes_simulate_cpp(Rcpp::NumericVector mu, Rcpp::NumericMatrix alpha,
    Rcpp::NumericMatrix beta, Rcpp::NumericMatrix initial, Rcpp::List mapping,
    Rcpp::NumericMatrix N, Rcpp::NumericMatrix Nc, Rcpp::NumericVector dt,
    Rcpp::NumericVector type, Rcpp::NumericVector mark) {
    const int d = beta.nrow(), k = beta.ncol(), size = dt.size();
    if (d < 1 || k < 1 || k > INT_MAX / d || size < 1 || mu.size() != d ||
        alpha.nrow() != d || alpha.ncol() != k || initial.nrow() != d ||
        initial.ncol() != k || N.nrow() != size || Nc.nrow() != size ||
        N.ncol() != d || Nc.ncol() != d || type.size() != size ||
        mark.size() != size || mapping.size() != d)
        Rcpp::stop("Invalid native simulation dimensions.");
    std::vector<std::vector<bool> > selected(d, std::vector<bool>(k, false));
    for (int i = 0; i < d; ++i) {
        if (!R_FINITE(mu[i]) || mu[i] <= 0.0)
            Rcpp::stop("Native simulation requires positive baselines.");
        Rcpp::IntegerVector columns = Rcpp::as<Rcpp::IntegerVector>(mapping[i]);
        for (int j = 0; j < columns.size(); ++j) {
            if (columns[j] < 0 || columns[j] >= k)
                Rcpp::stop("Invalid native simulation kernel mapping.");
            selected[i][columns[j]] = true;
        }
    }
    for (int j = 0; j < k; ++j) for (int i = 0; i < d; ++i) {
        if (!R_FINITE(alpha(i, j)) || !R_FINITE(beta(i, j)) || beta(i, j) <= 0.0 ||
            beta(i, j) != beta(i, 0))
            Rcpp::stop("Native simulation requires finite parameters and equal within-row decays.");
    }
    std::vector<double> state(initial.begin(), initial.end());
    if (!valid_state(state, d, k)) Rcpp::stop("Invalid native simulation initial state.");
    // R initially aliases type <- inter_arrival <- numeric(size). Its first
    // assignment separates them through copy-on-write; do the same explicitly.
    type = Rcpp::clone(type);
    Rcpp::NumericMatrix lambda(size, d), post(size, d), components(size, d * k), after(size, d * k);
    for (int i = 0; i < d; ++i) {
        lambda(0, i) = post(0, i) = mu[i] + row_sum(state, i, d, k);
        for (int j = 0; j < k; ++j) components(0, i * k + j) = after(0, i * k + j) = state[i + d * j];
    }
    int generated = 1;
    if (size > 1) {
        Rcpp::RNGScope scope;
        std::vector<double> candidates(d), uniforms(d), decayed(d * k);
        for (int n = 1; n < size; ++n) {
            if ((n & 1023) == 0) Rcpp::checkUserInterrupt();
            if (!valid_state(state, d, k))
                Rcpp::stop("Simulation requires finite, nonnegative baselines and nonnegative, decreasing excitation; signed components must satisfy the decay-rate cumulative-sum condition (see ?hsim).");
            // Preserve stats::rexp(D, rate=mu), then stats::runif(D): never interleave.
            for (int i = 0; i < d; ++i) candidates[i] = R::rexp(1.0 / mu[i]);
            for (int i = 0; i < d; ++i) uniforms[i] = R::runif(0.0, 1.0);
            double wait = R_PosInf;
            int event = -1;
            // The original cbind is searched in column-major order; ties retain
            // the first baseline row, then the first excitation row.
            for (int i = 0; i < d; ++i) if (candidates[i] < wait) {
                wait = candidates[i]; event = i;
            }
            for (int i = 0; i < d; ++i) {
                const double sum = row_sum(state, i, d, k);
                const double scaled = beta(i, 0) * std::log(uniforms[i]) / sum;
                if (sum > 0.0 && scaled > -1.0) {
                    const double excited = -std::log1p(scaled) / beta(i, 0);
                    if (excited < wait) { wait = excited; event = i; }
                }
            }
            if (!R_FINITE(wait)) {
                Rcpp::warning("Simulation stopped early: no further event occurs on this path.");
                break;
            }
            dt[n] = wait; type[n] = event + 1;
            for (int i = 0; i < d; ++i) {
                N(n, i) = N(n - 1, i); Nc(n, i) = Nc(n - 1, i);
            }
            N(n, event) += 1.0; Nc(n, event) += mark[n];
            for (int j = 0; j < k; ++j) for (int i = 0; i < d; ++i) {
                const int index = i + d * j;
                decayed[index] = state[index] * std::exp(-beta(i, j) * wait);
                components(n, i * k + j) = decayed[index];
            }
            for (int i = 0; i < d; ++i) lambda(n, i) = mu[i] + row_sum(decayed, i, d, k);
            for (int j = 0; j < k; ++j) for (int i = 0; i < d; ++i) {
                const int index = i + d * j;
                // Even unselected columns add zero, as in the R matrix addition.
                state[index] = decayed[index] + (selected[event][j] ? alpha(i, j) : 0.0);
                after(n, i * k + j) = state[index];
            }
            for (int i = 0; i < d; ++i) post(n, i) = mu[i] + row_sum(state, i, d, k);
            generated = n + 1;
        }
    }
    return Rcpp::List::create(Rcpp::Named("inter_arrival") = dt,
        Rcpp::Named("type") = type, Rcpp::Named("mark") = mark,
        Rcpp::Named("N") = N, Rcpp::Named("Nc") = Nc,
        Rcpp::Named("lambda") = lambda, Rcpp::Named("lambda_component") = components,
        Rcpp::Named("rambda") = post, Rcpp::Named("rambda_component") = after,
        Rcpp::Named("generated") = generated);
}
