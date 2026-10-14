#include <Rcpp.h>
#include <algorithm>
#include <climits>
#include <vector>
#include "exp_kernel.h"

// Numeric Hawkes kernel with optional precomputed observed-state impacts.
// No R callbacks, shared mutable state, or RNG calls.
// [[Rcpp::export(rng = false)]]
Rcpp::List hawkes_loglik_cpp(Rcpp::NumericVector dt, Rcpp::IntegerVector type,
    Rcpp::NumericVector mark, Rcpp::NumericVector mu, Rcpp::NumericMatrix alpha,
    Rcpp::NumericMatrix beta, Rcpp::Nullable<Rcpp::NumericMatrix> eta,
    Rcpp::NumericMatrix initial, Rcpp::List mapping, bool infer, bool gradient = false,
    Rcpp::Nullable<Rcpp::NumericMatrix> impact = R_NilValue) {

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
    Rcpp::NumericMatrix event_impact;
    const bool has_impact = impact.isNotNull();
    if (has_impact) {
        event_impact = Rcpp::NumericMatrix(impact.get());
        if (event_impact.nrow() != components || event_impact.ncol() != size)
            Rcpp::stop("Native impact dimensions must be components by observations.");
        for (R_xlen_t k = 0; k < event_impact.size(); ++k) {
            if ((k & 65535) == 0) Rcpp::checkUserInterrupt();
            if (!R_FINITE(event_impact[k])) Rcpp::stop("Invalid native impact value.");
        }
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
    // Local sensitivities keep score work O(events * components). R handles
    // parameter sharing and the chain rule for a parameter-dependent initial state.
    std::vector<double> da, db, de, dc;
    std::vector<long double> ga, gb, ge, gc, gm;
    if (gradient) {
        da.assign(components, 0.0); db.assign(components, 0.0);
        de.assign(components, 0.0); dc.assign(components, 1.0);
        ga.assign(components, 0.0L); gb.assign(components, 0.0L);
        ge.assign(components, 0.0L); gc.assign(components, 0.0L);
        gm.assign(dimensions, 0.0L);
    }
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
    std::vector<double> event_intensity;
    if (has_impact && gradient) event_intensity.resize(size);
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
                const double decay = std::exp(-scaled_time);
                if (gradient) {
                    const long double weight = scaled_time == 0.0 ? dt[n] :
                        -std::expm1(-scaled_time) / beta[index];
                    ga[index] -= da[index] * weight;
                    gb[index] -= db[index] * weight +
                        emhawkes::integrated_component_beta(state[index], beta[index], dt[n], decay, integral);
                    if (has_eta) ge[index] -= de[index] * weight;
                    gc[index] -= dc[index] * weight;
                    da[index] *= decay;
                    db[index] = db[index] * decay - dt[n] * (state[index] * decay);
                    if (has_eta) de[index] *= decay;
                    dc[index] *= decay;
                }
                state[index] *= decay;
                row_totals[i] += state[index];
                if (infer) lambda_component(n, i * kernels + j) = state[index];
            }
        }
        const int event_type = type[n] - 1;
        event_log += std::log(mu[event_type] + static_cast<double>(row_totals[event_type]));
        if (gradient) {
            const double intensity = mu[event_type] + static_cast<double>(row_totals[event_type]);
            const double denominator = intensity > 0.0 ? intensity : R_NaN;
            if (has_impact) event_intensity[n] = denominator;
            gm[event_type] += 1.0 / denominator;
            for (int j = 0; j < kernels; ++j) {
                const int index = event_type + dimensions * j;
                ga[index] += da[index] / denominator;
                gb[index] += db[index] / denominator;
                if (has_eta) ge[index] += de[index] / denominator;
                gc[index] += dc[index] / denominator;
            }
        }
        for (int i = 0; i < dimensions; ++i) {
            const double baseline = mu[i] * dt[n];
            if (gradient) gm[i] -= dt[n];
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
                if (has_impact) state[index] += event_impact(index, n);
                if (gradient) {
                    da[index] += 1.0;
                    if (has_eta) de[index] += mark[n] - 1.0;
                }
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
    Rcpp::List result = Rcpp::List::create(Rcpp::Named("loglikelihood") = event_log + integrated_log);
    if (gradient) {
        Rcpp::NumericVector score_mu(dimensions);
        Rcpp::NumericMatrix score_alpha(dimensions, kernels), score_beta(dimensions, kernels),
            score_eta(dimensions, kernels), score_initial(dimensions, kernels);
        const bool valid = R_FINITE(event_log + integrated_log);
        for (int i = 0; i < dimensions; ++i) score_mu[i] = valid ? static_cast<double>(gm[i]) : R_NaN;
        for (int j = 0; j < components; ++j) {
            score_alpha[j] = valid ? static_cast<double>(ga[j]) : R_NaN;
            score_beta[j] = valid ? static_cast<double>(gb[j]) : R_NaN;
            score_eta[j] = valid ? static_cast<double>(ge[j]) : R_NaN;
            score_initial[j] = valid ? static_cast<double>(gc[j]) : R_NaN;
        }
        result["score"] = Rcpp::List::create(Rcpp::Named("mu") = score_mu,
            Rcpp::Named("alpha") = score_alpha, Rcpp::Named("beta") = score_beta,
            Rcpp::Named("eta") = score_eta, Rcpp::Named("initial") = score_initial);
    }
    if (has_impact && gradient) {
        // Reverse sensitivity to each post-event jump. No events-by-parameters
        // tensor is needed; R contracts these scores with callback Jacobians.
        Rcpp::NumericMatrix score_impact(components, size);
        const bool valid = R_FINITE(event_log + integrated_log);
        if (!valid) std::fill(score_impact.begin(), score_impact.end(), R_NaN);
        else {
            std::vector<long double> adjoint(components, 0.0L);
            for (int n = size - 1; n >= 1; --n) {
                if ((n & 1023) == 0) Rcpp::checkUserInterrupt();
                const int event_type = type[n] - 1;
                const std::vector<int>& selected = columns[event_type];
                for (std::size_t j = 0; j < selected.size(); ++j)
                    for (int i = 0; i < dimensions; ++i) {
                        const int index = i + dimensions * selected[j];
                        score_impact(index, n) = static_cast<double>(adjoint[index]);
                    }
                for (int j = 0; j < kernels; ++j)
                    for (int i = 0; i < dimensions; ++i) {
                        const int index = i + dimensions * j;
                        const double x = beta[index] * dt[n];
                        const double decay = std::exp(-x);
                        const long double weight = x == 0.0 ? dt[n] : -std::expm1(-x) / beta[index];
                        adjoint[index] = adjoint[index] * decay - weight;
                        if (i == event_type) adjoint[index] += decay / event_intensity[n];
                    }
            }
        }
        result["impact_score"] = score_impact;
    }
    if (!infer) return result;
    Rcpp::List residual_output(dimensions);
    for (int i = 0; i < dimensions; ++i) residual_output[i] = Rcpp::wrap(residuals[i]);
    result["lambda"] = lambda;
    result["lambda_component"] = lambda_component;
    result["rambda"] = rambda;
    result["rambda_component"] = rambda_component;
    result["residuals"] = residual_output;
    return result;
}
