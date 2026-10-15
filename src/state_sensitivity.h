#ifndef EMHAWKES_STATE_SENSITIVITY_H
#define EMHAWKES_STATE_SENSITIVITY_H
#include <Rcpp.h>
#include <climits>
#include <memory>
#include <vector>
#include "jacobian_call.h"
#include "exp_kernel.h"

namespace emhawkes {
class StateSensitivity {
    int dimensions_, components_, parameters_;
    Rcpp::CharacterVector names_;
    Rcpp::NumericMatrix alpha_, beta_, eta_, mu_, baseline_;
    bool has_baseline_;
    std::vector<long double> state_, jump_, score_;
    std::unique_ptr<JacobianCall> impact_param_, impact_state_;
    Rcpp::NumericMatrix matrix(Rcpp::List input, const char* name, int rows) {
        SEXP value = input[name];
        if ((TYPEOF(value) != REALSXP && TYPEOF(value) != INTSXP) || !Rf_isMatrix(value) || Rf_isObject(value))
            Rcpp::stop("Invalid native sensitivity matrix: %s", name);
        Rcpp::NumericMatrix result = Rcpp::as<Rcpp::NumericMatrix>(value);
        if (result.nrow() != rows || result.ncol() != parameters_)
            Rcpp::stop("Invalid native sensitivity dimensions: %s", name);
        for (R_xlen_t k = 0; k < result.size(); ++k)
            if (!R_FINITE(result[k])) Rcpp::stop("Nonfinite native sensitivity: %s", name);
        return result;
    }
public:
    bool stateful;
    StateSensitivity(Rcpp::List input, int dimensions, int kernels, int size, bool baseline,
                     bool impact) : dimensions_(dimensions), components_(dimensions * kernels),
                     names_(input["names"]), has_baseline_(baseline) {
        parameters_ = names_.size();
        if (parameters_ && components_ > INT_MAX / parameters_)
            Rcpp::stop("Native sensitivity dimensions are too large.");
        mu_ = matrix(input, "mu", dimensions_); alpha_ = matrix(input, "alpha", components_);
        beta_ = matrix(input, "beta", components_); eta_ = matrix(input, "eta", components_);
        Rcpp::NumericMatrix initial = matrix(input, "initial", components_);
        state_.assign(initial.begin(), initial.end()); jump_.resize(state_.size());
        score_.assign(parameters_, 0.0L);
        stateful = Rcpp::as<bool>(input["stateful"]);
        if (has_baseline_) {
            SEXP b = input["baseline"];
            if ((TYPEOF(b) != REALSXP && TYPEOF(b) != INTSXP) || !Rf_isMatrix(b) || Rf_isObject(b))
                Rcpp::stop("Invalid native baseline Jacobian path.");
            baseline_ = Rcpp::as<Rcpp::NumericMatrix>(b);
            if (baseline_.nrow() != dimensions_ * parameters_ ||
                (baseline_.ncol() != 1 && baseline_.ncol() != size))
                Rcpp::stop("Invalid native baseline Jacobian path dimensions.");
            for (R_xlen_t k = 0; k < baseline_.size(); ++k)
                if (!R_FINITE(baseline_[k])) Rcpp::stop("Nonfinite native baseline Jacobian.");
        }
        if (impact) {
            Rcpp::List arguments = input["impact_arguments"];
            impact_param_.reset(new JacobianCall(input["impact_param"], arguments["param"],
                components_, parameters_, names_, true, stateful, "jacobian$impact$param"));
            if (stateful) impact_state_.reset(new JacobianCall(input["impact_state"], arguments["state"],
                components_, components_, Rcpp::CharacterVector(0), false, true, "jacobian$impact$state"));
        } else if (stateful) Rcpp::stop("State feedback requires an impact callback.");
    }
    double mu(int row, int parameter, int n) const {
        return has_baseline_ ? baseline_(row + dimensions_ * parameter,
            baseline_.ncol() == 1 ? 0 : n) : mu_(row, parameter);
    }
    void decay(int index, double old_state, double rate, double dt, double scaled,
               double decay, double integral) {
        const long double weight = scaled == 0.0 ? dt : -std::expm1(-scaled) / rate;
        const long double beta_partial = integrated_component_beta(old_state, rate, dt, decay, integral);
        for (int p = 0; p < parameters_; ++p) {
            const int k = index + components_ * p;
            const long double before = state_[k];
            const double db = beta_(index, p);
            if (before != 0.0L) score_[p] -= before * weight;
            if (db != 0.0) score_[p] -= beta_partial * db;
            state_[k] = before * decay;
            if (db != 0.0) state_[k] -= (static_cast<long double>(dt) * (old_state * decay)) * db;
        }
    }
    void event(int n, int event_type, double intensity, double dt) {
        const double denominator = intensity > 0.0 ? intensity : R_NaN;
        for (int p = 0; p < parameters_; ++p) {
            long double derivative = mu(event_type, p, n);
            for (int index = event_type; index < components_; index += dimensions_)
                derivative += state_[index + components_ * p];
            score_[p] += derivative / denominator;
            for (int i = 0; i < dimensions_; ++i) {
                const double dmu = mu(i, p, n);
                if (dmu != 0.0) score_[p] -= static_cast<long double>(dt) * dmu;
            }
        }
    }
    void jump(int n, double mark, const std::vector<int>& selected,
              const Rcpp::NumericMatrix& beta, const std::vector<double>& state) {
        if (impact_param_) impact_param_->evaluate(n + 1, beta, state);
        if (impact_state_) impact_state_->evaluate(n + 1, beta, state);
        std::fill(jump_.begin(), jump_.end(), 0.0L);
        for (std::size_t j = 0; j < selected.size(); ++j)
            for (int i = 0; i < dimensions_; ++i) {
                const int index = i + dimensions_ * selected[j];
                for (int p = 0; p < parameters_; ++p) {
                    long double value = alpha_(index, p);
                    if (eta_(index, p) != 0.0) value += eta_(index, p) * (mark - 1.0);
                    if (impact_param_) value += impact_param_->at(index, p);
                    jump_[index + components_ * p] = value;
                }
                if (impact_state_) for (int k = 0; k < components_; ++k) {
                    const double weight = impact_state_->at(index, k);
                    if (weight == 0.0) continue;
                    for (int p = 0; p < parameters_; ++p)
                        jump_[index + components_ * p] += weight * state_[k + components_ * p];
                }
            }
        for (std::size_t k = 0; k < state_.size(); ++k) state_[k] += jump_[k];
    }
    Rcpp::NumericVector gradient(bool valid) const {
        Rcpp::NumericVector result(parameters_);
        for (int p = 0; p < parameters_; ++p) result[p] = valid ? static_cast<double>(score_[p]) : R_NaN;
        result.names() = names_; return result;
    }
};
}
#endif
