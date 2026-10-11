#include <Rcpp.h>
#include <algorithm>
#include <climits>
#include <vector>
#include "exp_kernel.h"

// Numeric Hawkes only: no R callbacks, shared mutable state, or RNG calls.
// [[Rcpp::export(rng = false)]]
Rcpp::List hawkes_loglik_cpp(Rcpp::NumericVector dt, Rcpp::IntegerVector type,
    Rcpp::NumericVector mark, Rcpp::NumericVector mu, Rcpp::NumericMatrix alpha,
    Rcpp::NumericMatrix beta, Rcpp::Nullable<Rcpp::NumericMatrix> eta,
    Rcpp::NumericMatrix initial, Rcpp::List mapping, bool infer) {

    const int dimensions = beta.nrow(), kernels = beta.ncol();
    if (dimensions < 1 || kernels < 1 || dt.size() < 1 || dt.size() > INT_MAX ||
        kernels > INT_MAX / dimensions) Rcpp::stop("Invalid native Hawkes dimensions.");
    const int size = static_cast<int>(dt.size()), components = dimensions * kernels;
    if (type.size() != size || mark.size() != size || mu.size() != dimensions ||
        alpha.nrow() != dimensions || alpha.ncol() != kernels ||
        initial.nrow() != dimensions || initial.ncol() != kernels ||
        mapping.size() != dimensions) Rcpp::stop("Native Hawkes inputs have inconsistent dimensions.");
    Rcpp::NumericMatrix mark_effect;
    const bool has_eta = eta.isNotNull();
    if (has_eta) {
        mark_effect = Rcpp::NumericMatrix(eta.get());
        if (mark_effect.nrow() != dimensions || mark_effect.ncol() != kernels)
            Rcpp::stop("Native eta dimensions do not match beta.");
    }
    for (int n = 0; n < size; ++n) {
        if ((n & 65535) == 0) Rcpp::checkUserInterrupt();
        if (!R_FINITE(dt[n]) || dt[n] < 0 || !R_FINITE(mark[n]) ||
            (n > 0 && (type[n] < 1 || type[n] > dimensions)))
            Rcpp::stop("Invalid native Hawkes observation.");
    }
    if (dt[0] != 0) Rcpp::stop("Native inter-arrival times must start with zero.");
    for (int i = 0; i < dimensions; ++i) {
        if (!R_FINITE(mu[i]) || mu[i] < 0) Rcpp::stop("Invalid native baseline intensity.");
    }
    for (int j = 0; j < components; ++j) {
        if ((j & 65535) == 0) Rcpp::checkUserInterrupt();
        if (!R_FINITE(beta[j]) || beta[j] <= 0 || !R_FINITE(alpha[j]) ||
            !R_FINITE(initial[j]) || (has_eta && !R_FINITE(mark_effect[j])))
            Rcpp::stop("Invalid native Hawkes parameter or initial state.");
    }
    std::vector<std::vector<int> > columns(dimensions);
    for (int i = 0; i < dimensions; ++i) {
        Rcpp::IntegerVector selected = Rcpp::as<Rcpp::IntegerVector>(mapping[i]);
        if (selected.size() < 1) Rcpp::stop("Empty native kernel mapping.");
        std::vector<bool> seen(kernels, false);
        for (int j = 0; j < selected.size(); ++j) {
            const int column = selected[j];
            if (column < 0 || column >= kernels || seen[column])
                Rcpp::stop("Invalid native kernel mapping.");
            seen[column] = true;
            columns[i].push_back(column);
        }
    }

    // R inputs are read only. Recursive state owns its storage.
    std::vector<double> state(initial.begin(), initial.end());
    std::vector<long double> row_integrals(dimensions), row_totals(dimensions);
    std::vector<long double> residual_integrals(dimensions, 0.0L);
    std::vector<bool> observed(dimensions, false);
    std::vector<std::vector<double> > residuals(dimensions);
    Rcpp::NumericMatrix lambda, lambda_component, rambda, rambda_component;
    if (infer) {
        lambda = Rcpp::NumericMatrix(size, dimensions);
        rambda = Rcpp::NumericMatrix(size, dimensions);
        lambda_component = Rcpp::NumericMatrix(size, components);
        rambda_component = Rcpp::NumericMatrix(size, components);
        for (int i = 0; i < dimensions; ++i) {
            long double total = 0.0L;
            for (int j = 0; j < kernels; ++j) {
                const double value = state[i + dimensions * j];
                total += value;
                const int output_column = i * kernels + j;
                lambda_component(0, output_column) = value;
                rambda_component(0, output_column) = value;
            }
            lambda(0, i) = rambda(0, i) = mu[i] + static_cast<double>(total);
        }
    }
    double event_log = 0.0, integrated_log = 0.0;
    for (int n = 1; n < size; ++n) {
        if ((n & 1023) == 0) Rcpp::checkUserInterrupt();
        std::fill(row_integrals.begin(), row_integrals.end(), 0.0L);
        std::fill(row_totals.begin(), row_totals.end(), 0.0L);
        long double total_integral = 0.0L, baseline_integral = 0.0L;
        for (int j = 0; j < kernels; ++j) {
            for (int i = 0; i < dimensions; ++i) {
                const int index = i + dimensions * j;
                const double scaled_time = beta[index] * dt[n];
                const double integral = emhawkes::integrated_component(state[index], beta[index], dt[n], scaled_time);
                row_integrals[i] += integral;
                total_integral += integral;
                state[index] *= std::exp(-scaled_time);
                row_totals[i] += state[index];
                if (infer) lambda_component(n, i * kernels + j) = state[index];
            }
        }
        const int event_type = type[n] - 1;
        event_log += std::log(mu[event_type] + static_cast<double>(row_totals[event_type]));
        for (int i = 0; i < dimensions; ++i) {
            const double baseline = mu[i] * dt[n];
            baseline_integral += baseline;
            if (infer) {
                const double phi = static_cast<double>(row_integrals[i]) + baseline;
                integrated_log -= phi;
                residual_integrals[i] += phi;
                lambda(n, i) = mu[i] + static_cast<double>(row_totals[i]);
            }
        }
        if (!infer) {
            // Preserve R's separate sums and per-event accumulation order.
            integrated_log = integrated_log - static_cast<double>(total_integral) - static_cast<double>(baseline_integral);
        }
        if (infer) {
            if (observed[event_type]) residuals[event_type].push_back(static_cast<double>(residual_integrals[event_type]));
            observed[event_type] = true;
            residual_integrals[event_type] = 0.0L;
        }
        const std::vector<int>& selected = columns[event_type];
        for (std::size_t j = 0; j < selected.size(); ++j) {
            for (int i = 0; i < dimensions; ++i) {
                const int index = i + dimensions * selected[j];
                state[index] += alpha[index];
                if (has_eta) state[index] += mark_effect[index] * (mark[n] - 1.0);
            }
        }
        if (infer) {
            for (int i = 0; i < dimensions; ++i) {
                long double total = 0.0L;
                for (int j = 0; j < kernels; ++j) {
                    const double value = state[i + dimensions * j];
                    total += value;
                    rambda_component(n, i * kernels + j) = value;
                }
                rambda(n, i) = mu[i] + static_cast<double>(total);
            }
        }
    }
    if (!infer) return Rcpp::List::create(Rcpp::Named("loglikelihood") = event_log + integrated_log);
    Rcpp::List residual_output(dimensions);
    for (int i = 0; i < dimensions; ++i) residual_output[i] = Rcpp::wrap(residuals[i]);
    return Rcpp::List::create(Rcpp::Named("loglikelihood") = event_log + integrated_log,
        Rcpp::Named("lambda") = lambda, Rcpp::Named("lambda_component") = lambda_component,
        Rcpp::Named("rambda") = rambda, Rcpp::Named("rambda_component") = rambda_component,
        Rcpp::Named("residuals") = residual_output);
}
