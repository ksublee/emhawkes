#include <Rcpp.h>
#include <climits>

// Evaluate an observed, piecewise-constant baseline without an R event loop.
// Every call owns its argument cells; retained promises keep their event values.
// [[Rcpp::export(rng = false)]]
Rcpp::List hawkes_baseline_cpp(Rcpp::Function callback, Rcpp::List arguments,
    int dimensions, int size, Rcpp::Nullable<Rcpp::NumericMatrix> initial = R_NilValue,
    bool infer = false, bool integer_n = true) {
    if (dimensions < 1 || size < 1 || size == INT_MAX)
        Rcpp::stop("Invalid native baseline path dimensions.");
    const bool has_initial = initial.isNotNull();
    const int columns = size + (infer && size > 1 ? 1 : 0);
    Rcpp::NumericMatrix path(dimensions, columns);
    if (has_initial) {
        Rcpp::NumericMatrix first(initial.get());
        if (first.nrow() != dimensions || first.ncol() != 1)
            Rcpp::stop("Invalid native initial baseline dimensions.");
        for (int i = 0; i < dimensions; ++i) path(i, 0) = first[i];
    }
    if (arguments.size() && Rf_isNull(arguments.names())) Rcpp::stop("Baseline arguments must be named.");
    Rcpp::CharacterVector names = arguments.size() ? arguments.names() : Rcpp::CharacterVector(0);
    if (names.size() != arguments.size()) Rcpp::stop("Invalid baseline argument names.");
    Rcpp::RObject cells = Rf_allocList(arguments.size());
    Rcpp::RObject call = Rf_lcons(callback, cells);
    SEXP node = CDR(call), n_node = R_NilValue;
    for (int k = 0; k < arguments.size(); ++k, node = CDR(node)) {
        if (names[k] == NA_STRING) Rcpp::stop("Invalid baseline argument name.");
        const std::string name = Rcpp::as<std::string>(names[k]);
        if (name.empty()) Rcpp::stop("Invalid baseline argument name.");
        SET_TAG(node, Rf_install(name.c_str()));
        SETCAR(node, arguments[k]);
        if (name == "n") n_node = node;
    }
    // Capability expressions in unavailable latent arguments belong to this
    // namespace; callback closures retain their own lexical environments.
    Rcpp::Environment environment = Rcpp::Environment::namespace_env("emhawkes");
    const int last = size + (size > 1 ? 1 : 0);
    for (int n = has_initial ? 2 : 1; n <= last; ++n) {
        if ((n & 1023) == 0) Rcpp::checkUserInterrupt();
        if (n_node != R_NilValue)
            SETCAR(n_node, integer_n ? Rf_ScalarInteger(n) : Rf_ScalarReal(n));
        Rcpp::RObject value = Rcpp::Rcpp_fast_eval(call, environment);
        // R evaluates the terminal next baseline even when inference is off.
        // Its unused result has no shape or finite-value requirement.
        if (n > columns) continue;
        if ((TYPEOF(value) != REALSXP && TYPEOF(value) != INTSXP) ||
            !Rf_isMatrix(value) || Rf_isObject(value))
            return Rcpp::List::create(Rcpp::Named("unsupported") =
                "Observed mu must return a plain numeric dimensions-by-one matrix.");
        Rcpp::NumericMatrix mu = Rcpp::as<Rcpp::NumericMatrix>(value);
        if (mu.nrow() != dimensions || mu.ncol() != 1)
            return Rcpp::List::create(Rcpp::Named("unsupported") =
                "Observed mu must return a plain numeric dimensions-by-one matrix.");
        for (int i = 0; i < dimensions; ++i) {
            const double x = mu[i];
            if (!R_FINITE(x) || x < 0)
                return Rcpp::List::create(Rcpp::Named("invalid") = true);
            path(i, n - 1) = x;
        }
    }
    return Rcpp::List::create(Rcpp::Named("path") = path);
}
