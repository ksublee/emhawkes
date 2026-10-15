#ifndef EMHAWKES_JACOBIAN_CALL_H
#define EMHAWKES_JACOBIAN_CALL_H
#include <Rcpp.h>
#include <algorithm>
#include <string>
#include <vector>

namespace emhawkes {
// Own the call cells and each returned matrix. Changing cells never changes
// argument values retained by a callback. A constant Jacobian is validated once.
class JacobianCall {
    Rcpp::RObject call_;
    Rcpp::NumericMatrix value_;
    Rcpp::CharacterVector parameter_names_;
    std::vector<int> order_;
    SEXP n_node_, current_node_;
    int rows_, cols_;
    bool function_, named_, current_available_;
    std::string label_;
    void validate(SEXP value) {
        if ((TYPEOF(value) != REALSXP && TYPEOF(value) != INTSXP) ||
            !Rf_isMatrix(value) || Rf_isObject(value))
            Rcpp::stop(label_ + " must return a finite plain numeric matrix.");
        value_ = Rcpp::as<Rcpp::NumericMatrix>(value);
        if (value_.nrow() != rows_ || value_.ncol() != cols_)
            Rcpp::stop(label_ + " has incorrect dimensions.");
        for (R_xlen_t k = 0; k < value_.size(); ++k) {
            if ((k & 65535) == 0) Rcpp::checkUserInterrupt();
            if (!R_FINITE(value_[k])) Rcpp::stop(label_ + " must contain finite values.");
        }
        order_.resize(cols_);
        for (int p = 0; p < cols_; ++p) order_[p] = p;
        if (named_ && cols_) {
            SEXP dn = Rf_getAttrib(value_, R_DimNamesSymbol);
            if (Rf_isNull(dn) || Rf_isNull(VECTOR_ELT(dn, 1)))
                Rcpp::stop(label_ + " column names must match all fitted parameters.");
            Rcpp::CharacterVector names = VECTOR_ELT(dn, 1);
            if (names.size() != cols_) Rcpp::stop(label_ + " has invalid column names.");
            std::vector<bool> used(cols_, false);
            for (int p = 0; p < cols_; ++p) {
                int found = -1;
                for (int j = 0; j < cols_; ++j)
                    if (names[j] != NA_STRING && names[j] == parameter_names_[p]) {
                        if (found >= 0) Rcpp::stop(label_ + " has duplicate column names.");
                        found = j;
                    }
                if (found < 0 || used[found])
                    Rcpp::stop(label_ + " column names must match all fitted parameters.");
                used[found] = true; order_[p] = found;
            }
        }
    }
public:
    JacobianCall(SEXP value, Rcpp::List arguments, int rows, int cols,
                 Rcpp::CharacterVector names, bool named, bool current_available,
                 const std::string& label) : parameter_names_(names), n_node_(R_NilValue),
                 current_node_(R_NilValue), rows_(rows), cols_(cols),
                 function_(Rf_isFunction(value)), named_(named),
                 current_available_(current_available), label_(label) {
        if (!function_) { validate(value); return; }
        if (arguments.size() && Rf_isNull(arguments.names()))
            Rcpp::stop("Jacobian arguments must be named.");
        Rcpp::CharacterVector argument_names = arguments.size() ? arguments.names() : Rcpp::CharacterVector(0);
        if (argument_names.size() != arguments.size()) Rcpp::stop("Invalid Jacobian argument names.");
        Rcpp::RObject cells = Rf_allocList(arguments.size());
        call_ = Rf_lcons(value, cells);
        SEXP node = CDR(call_);
        for (int k = 0; k < arguments.size(); ++k, node = CDR(node)) {
            if (argument_names[k] == NA_STRING) Rcpp::stop("Invalid Jacobian argument name.");
            const std::string name = Rcpp::as<std::string>(argument_names[k]);
            if (name.empty()) Rcpp::stop("Invalid Jacobian argument name.");
            SET_TAG(node, Rf_install(name.c_str())); SETCAR(node, arguments[k]);
            if (name == "n") n_node_ = node;
            if (name == "lambda_component_n") current_node_ = node;
        }
    }
    void evaluate(int n, const Rcpp::NumericMatrix& beta, const std::vector<double>& state) {
        if (!function_) return;
        if (n_node_ != R_NilValue) SETCAR(n_node_, Rf_ScalarInteger(n));
        if (current_node_ != R_NilValue && current_available_) {
            Rcpp::NumericMatrix current = Rcpp::clone(beta);
            std::copy(state.begin(), state.end(), current.begin()); SETCAR(current_node_, current);
        }
        Rcpp::RObject result = Rcpp::Rcpp_fast_eval(call_, Rcpp::Environment::namespace_env("emhawkes"));
        if (n_node_ != R_NilValue) SETCAR(n_node_, R_NilValue);
        if (current_node_ != R_NilValue && current_available_) SETCAR(current_node_, R_NilValue);
        validate(result);
    }
    double at(int row, int column) const { return value_(row, order_[column]); }
};
}
#endif
