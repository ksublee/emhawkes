#include <Rcpp.h>
#include <algorithm>
#include <climits>
#include <string>
#include <vector>
#include "exp_kernel.h"
#include "state_sensitivity.h"

// Numeric Hawkes kernel with optional precomputed observed-state impacts.
// Optional observed baselines and sequential impacts support numerical replay.
// Recursive storage is private; this kernel never draws random numbers.
// [[Rcpp::export(rng = false)]]
Rcpp::List hawkes_loglik_cpp(Rcpp::NumericVector dt, Rcpp::IntegerVector type,
    Rcpp::NumericVector mark, Rcpp::NumericVector mu, Rcpp::NumericMatrix alpha,
    Rcpp::NumericMatrix beta, Rcpp::Nullable<Rcpp::NumericMatrix> eta,
    Rcpp::NumericMatrix initial, Rcpp::List mapping, bool infer, bool gradient = false,
    Rcpp::Nullable<Rcpp::NumericMatrix> impact = R_NilValue,
    Rcpp::Nullable<Rcpp::Function> impact_callback = R_NilValue,
    bool keep_lambda = false, bool keep_components = false,
    Rcpp::Nullable<Rcpp::NumericMatrix> baseline = R_NilValue,
    Rcpp::Nullable<Rcpp::List> impact_arguments = R_NilValue,
    Rcpp::Nullable<Rcpp::List> sensitivity = R_NilValue, bool strict_callback = false) {

    const bool full_gradient = sensitivity.isNotNull();
    const bool local_gradient = gradient && !full_gradient;
    const bool strict = full_gradient || strict_callback;
    if (full_gradient && (!gradient || infer))
        Rcpp::stop("Full state sensitivities require gradient = true and infer = false.");
    if (full_gradient && (keep_lambda || keep_components))
        Rcpp::stop("Analytic state sensitivities do not support full intensity histories.");

    const bool has_callback = impact_callback.isNotNull();
    if (has_callback != impact_arguments.isNotNull())
        Rcpp::stop("Sequential impacts require a callback and named arguments together.");
    if (has_callback && (local_gradient || impact.isNotNull()))
        Rcpp::stop("Sequential impacts do not support gradients or precomputed impacts.");
    if (!has_callback && (keep_lambda || keep_components))
        Rcpp::stop("State histories require a sequential impact callback.");

    const int dimensions = beta.nrow(), kernels = beta.ncol();
    if (dimensions < 1 || kernels < 1 || dt.size() < 1 || dt.size() > INT_MAX ||
        kernels > INT_MAX / dimensions) Rcpp::stop("Invalid native Hawkes dimensions.");
    const int size = static_cast<int>(dt.size()), components = dimensions * kernels;
    if (type.size() != size || mark.size() != size || mu.size() != dimensions ||
        alpha.nrow() != dimensions || alpha.ncol() != kernels ||
        initial.nrow() != dimensions || initial.ncol() != kernels ||
        mapping.size() != dimensions) Rcpp::stop("Native Hawkes inputs have inconsistent dimensions.");
    const bool has_baseline = baseline.isNotNull();
    Rcpp::NumericMatrix baseline_values;
    if (has_baseline) {
        if (local_gradient) Rcpp::stop("Observed baseline paths do not support gradients.");
        baseline_values = Rcpp::NumericMatrix(baseline.get());
        const R_xlen_t expected_columns = static_cast<R_xlen_t>(size) + (infer && size > 1 ? 1 : 0);
        if (baseline_values.nrow() != dimensions || baseline_values.ncol() != expected_columns)
            Rcpp::stop("Native baseline dimensions must include the terminal next baseline for inference.");
        for (R_xlen_t k = 0; k < baseline_values.size(); ++k) {
            if ((k & 65535) == 0) Rcpp::checkUserInterrupt();
            if (!R_FINITE(baseline_values[k]) || baseline_values[k] < 0)
                Rcpp::stop("Invalid native observed baseline.");
        }
    }
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

    std::unique_ptr<emhawkes::StateSensitivity> full;
    if (full_gradient) full.reset(new emhawkes::StateSensitivity(Rcpp::List(sensitivity.get()),
        dimensions, kernels, size, has_baseline, has_callback || has_impact));
    bool expose_current = !full_gradient || full->stateful;

    // R inputs are read only. Recursive state owns its storage.
    std::vector<double> state(initial.begin(), initial.end());
    // Local sensitivities keep score work O(events * components). R handles
    // parameter sharing and the chain rule for a parameter-dependent initial state.
    std::vector<double> da, db, de, dc;
    std::vector<long double> ga, gb, ge, gc, gm;
    if (local_gradient) {
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
    if (has_callback) {
        if (keep_lambda) {
            if (!infer) lambda = Rcpp::NumericMatrix(size, dimensions);
            Rcpp::CharacterVector names(dimensions);
            for (int i = 0; i < dimensions; ++i) {
                long double total = 0.0L;
                for (int j = 0; j < kernels; ++j) total += state[i + dimensions * j];
                lambda(0, i) = mu[i] + static_cast<double>(total);
                names[i] = "lambda" + std::to_string(i + 1);
            }
            Rcpp::colnames(lambda) = names;
        }
        if (keep_components) {
            if (!infer) lambda_component = Rcpp::NumericMatrix(size, components);
            Rcpp::CharacterVector names(components);
            for (int i = 0; i < dimensions; ++i)
                for (int j = 0; j < kernels; ++j) {
                    const int column = i * kernels + j;
                    lambda_component(0, column) = state[i + dimensions * j];
                    names[column] = "lambda" + std::to_string(i + 1) + std::to_string(j + 1);
                }
            Rcpp::colnames(lambda_component) = names;
        }
    }
    std::vector<double> event_intensity;
    if (has_impact && local_gradient) event_intensity.resize(size);
    // Build one named R call per proposal. Updating only the event/state nodes
    // avoids an R wrapper and do.call in every event. The call owns all nodes.
    Rcpp::RObject callback_call;
    SEXP n_node = R_NilValue, current_node = R_NilValue;
    SEXP lambda_node = R_NilValue, components_node = R_NilValue;
    if (has_callback) {
        Rcpp::List arguments(impact_arguments.get());
        if (arguments.size() && Rf_isNull(arguments.names())) Rcpp::stop("Sequential impact arguments must be named.");
        Rcpp::CharacterVector names = arguments.size() ? arguments.names() : Rcpp::CharacterVector(0);
        if (names.size() != arguments.size()) Rcpp::stop("Invalid sequential impact argument names.");
        Rcpp::RObject cells = Rf_allocList(arguments.size());
        callback_call = Rf_lcons(impact_callback.get(), cells);
        SEXP node = CDR(callback_call);
        for (int k = 0; k < arguments.size(); ++k, node = CDR(node)) {
            if (names[k] == NA_STRING) Rcpp::stop("Invalid sequential impact argument name.");
            const std::string name = Rcpp::as<std::string>(names[k]);
            if (name.empty()) Rcpp::stop("Invalid sequential impact argument name.");
            SET_TAG(node, Rf_install(name.c_str()));
            SETCAR(node, arguments[k]);
            if (name == "n") n_node = node;
            else if (name == "lambda_component_n") current_node = node;
            else if (name == "lambda") lambda_node = node;
            else if (name == "lambda_component") components_node = node;
        }
    }
    Rcpp::NumericMatrix callback_impact;
    // Analytic value-only passes leave unavailable current states as promises
    // that raise an error. Explicit current-state callbacks receive NULL here.
    if (strict_callback && current_node != R_NilValue && TYPEOF(CAR(current_node)) == LANGSXP)
        expose_current = false;
    double event_log = 0.0, integrated_log = 0.0;
    for (int n = 1; n < size; ++n) {
        if ((n & 1023) == 0) Rcpp::checkUserInterrupt();
        // R callbacks may retain histories. Copy before later writes, as R's
        // own indexed assignment would, so retained snapshots do not change.
        if (has_callback && keep_lambda && MAYBE_SHARED(lambda)) lambda = Rcpp::clone(lambda);
        if (has_callback && keep_components && MAYBE_SHARED(lambda_component))
            lambda_component = Rcpp::clone(lambda_component);
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
                if (local_gradient) {
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
                if (full_gradient) full->decay(index, state[index], beta[index], dt[n],
                    scaled_time, decay, integral);
                state[index] *= decay;
                row_totals[i] += state[index];
                if (infer || keep_components) lambda_component(n, i * kernels + j) = state[index];
            }
        }
        const int event_type = type[n] - 1;
        event_log += std::log((has_baseline ? baseline_values(event_type, n) : mu[event_type]) +
                              static_cast<double>(row_totals[event_type]));
        if (full_gradient) full->event(n, event_type,
            (has_baseline ? baseline_values(event_type, n) : mu[event_type]) +
            static_cast<double>(row_totals[event_type]), dt[n]);
        if (local_gradient) {
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
            const double mu_at_event = has_baseline ? baseline_values(i, n) : mu[i];
            const double baseline = mu_at_event * dt[n];
            if (local_gradient) gm[i] -= dt[n];
            baseline_integral += baseline;
            if (infer) {
                const double phi = static_cast<double>(row_integrals[i]) + baseline;
                integrated_log -= phi;
                residual_integrals[i] += phi;
            }
            if (infer || keep_lambda) lambda(n, i) = mu_at_event + static_cast<double>(row_totals[i]);
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
        if (has_callback) {
            if (n_node != R_NilValue) SETCAR(n_node, Rf_ScalarInteger(n + 1));
            if (current_node != R_NilValue && expose_current) {
                // R's elementwise decay carries beta's attributes to the current
                // component matrix, including dimnames used by callbacks.
                Rcpp::NumericMatrix current = Rcpp::clone(beta);
                std::copy(state.begin(), state.end(), current.begin());
                SETCAR(current_node, current);
            }
            if (lambda_node != R_NilValue && !strict)
                SETCAR(lambda_node, keep_lambda ? static_cast<SEXP>(lambda) : R_NilValue);
            if (components_node != R_NilValue && !strict)
                SETCAR(components_node, keep_components ? static_cast<SEXP>(lambda_component) : R_NilValue);
            Rcpp::RObject value = Rcpp::Rcpp_fast_eval(callback_call, strict ?
                static_cast<SEXP>(Rcpp::Environment::namespace_env("emhawkes")) : R_GlobalEnv);
            // Release our references before the next copy-on-write check.
            if (n_node != R_NilValue) SETCAR(n_node, R_NilValue);
            if (current_node != R_NilValue && expose_current) SETCAR(current_node, R_NilValue);
            if (lambda_node != R_NilValue && !strict) SETCAR(lambda_node, R_NilValue);
            if (components_node != R_NilValue && !strict) SETCAR(components_node, R_NilValue);
            if ((TYPEOF(value) != REALSXP && TYPEOF(value) != INTSXP) ||
                !Rf_isMatrix(value) || Rf_isObject(value))
                return Rcpp::List::create(Rcpp::Named("unsupported") =
                    "impact must return a plain numeric matrix matching the kernel dimensions.");
            callback_impact = Rcpp::as<Rcpp::NumericMatrix>(value);
            if (callback_impact.nrow() != dimensions || callback_impact.ncol() != kernels)
                return Rcpp::List::create(Rcpp::Named("unsupported") =
                    "impact must return a plain numeric matrix matching the kernel dimensions.");
            for (int k = 0; k < components; ++k)
                if (!R_FINITE(callback_impact[k]))
                    return Rcpp::List::create(Rcpp::Named("loglikelihood") = R_NegInf,
                                              Rcpp::Named("nonfinite_impact") = true);
        }
        const std::vector<int>& selected = columns[event_type];
        if (full_gradient) full->jump(n, mark[n], selected, beta, state);
        for (std::size_t j = 0; j < selected.size(); ++j) {
            for (int i = 0; i < dimensions; ++i) {
                const int index = i + dimensions * selected[j];
                state[index] += alpha[index];
                if (has_eta) state[index] += mark_effect[index] * (mark[n] - 1.0);
                if (has_impact) state[index] += event_impact(index, n);
                if (has_callback) state[index] += callback_impact[index];
                if (local_gradient) {
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
                const double post_mu = has_baseline ? baseline_values(i, n + 1) : mu[i];
                rambda(n, i) = post_mu + static_cast<double>(total);
            }
        }
    }
    Rcpp::List result = Rcpp::List::create(Rcpp::Named("loglikelihood") = event_log + integrated_log);
    if (full_gradient) result["gradient"] = full->gradient(R_FINITE(event_log + integrated_log));
    if (local_gradient) {
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
    if (has_impact && local_gradient) {
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
