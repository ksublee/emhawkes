#include "simulation_state.h"

// Compile the equal-rate loop separately so its optimization cannot change
// generated code for the existing general numeric loop.

// Fresh output/count buffers are prepared privately by simulate_hawkes().
// No user input matrix is mutated. Support is decided in R before this call.
// The explicit scope below avoids even creating .Random.seed for size=1.
// [[Rcpp::export(rng = false)]]
Rcpp::List hawkes_simulate_equal_cpp(Rcpp::NumericVector mu, Rcpp::NumericMatrix alpha,
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
    bool equal = true;
    for (int j = 0; j < k; ++j) for (int i = 0; i < d; ++i) {
        if (!R_FINITE(alpha(i, j)) || !R_FINITE(beta(i, j)) || beta(i, j) <= 0.0)
            Rcpp::stop("Native simulation requires finite parameters and positive decays.");
        equal = equal && beta(i, j) == beta(i, 0);
    }
    if (!equal) Rcpp::stop("Native equal-decay simulation requires equal within-row decays.");
    RateGroups groups(d);
    for (int i = 0; i < d; ++i) {
        std::vector<double> rates;
        for (int j = 0; j < k; ++j) rates.push_back(beta(i, j));
        std::sort(rates.begin(), rates.end());
        rates.erase(std::unique(rates.begin(), rates.end()), rates.end());
        groups[i].resize(rates.size());
        // Preserve original column order inside each rate group.
        for (int j = 0; j < k; ++j) {
            const int group = std::lower_bound(rates.begin(), rates.end(), beta(i, j)) - rates.begin();
            groups[i][group].push_back(j);
        }
    }
    std::vector<double> state(initial.begin(), initial.end());
    if (!valid_state(state, d, k, beta, groups)) Rcpp::stop("Invalid native simulation initial state.");
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
        SimulationRNGScope scope;
        std::vector<double> candidates(d), uniforms(d), decayed(d * k);
        for (int n = 1; n < size; ++n) {
            if ((n & 1023) == 0) Rcpp::checkUserInterrupt();
            if (!valid_state(state, d, k, beta, groups))
                Rcpp::stop("Simulation requires finite, nonnegative baselines and nonnegative, decreasing excitation; signed components must satisfy the decay-rate cumulative-sum condition (see ?hsim).");
            double wait = R_PosInf;
            int event = -1;
            {
                // Preserve stats::rexp(D, rate=mu), then stats::runif(D): never interleave.
                for (int i = 0; i < d; ++i) candidates[i] = R::rexp(1.0 / mu[i]);
                for (int i = 0; i < d; ++i) uniforms[i] = R::runif(0.0, 1.0);
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
            }
            if (!R_FINITE(wait)) {
                simulation_warning("Simulation stopped early: no further event occurs on this path. Returned " +
                    std::to_string(n) + " of " + std::to_string(size) + " requested rows (including the initial row).");
                break;
            }
            dt[n] = wait; type[n] = event + 1;
            for (int i = 0; i < d; ++i) {
                N(n, i) = N(n - 1, i); Nc(n, i) = Nc(n - 1, i);
            }
            N(n, event) += 1.0; Nc(n, event) += mark[n];
            // Arrival uniforms are no longer needed; reuse their buffer for
            // row decay factors, retaining component products and sums.
            for (int i = 0; i < d; ++i)
                uniforms[i] = std::exp(-beta(i, 0) * wait);
            for (int j = 0; j < k; ++j) for (int i = 0; i < d; ++i) {
                const int index = i + d * j;
                decayed[index] = state[index] * uniforms[i];
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
