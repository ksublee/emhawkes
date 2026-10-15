library(emhawkes)

evaluate <- function(engine, object, dt, type = NULL, mark = NULL, initial = .1,
                     infer = FALSE, ...) {
  previous <- options(emhawkes.engine = engine)
  on.exit(options(previous))
  suppressMessages(logLik(object, dt, type = type, mark = mark,
                         lambda_component0 = initial, infer = infer, ...))
}
close <- function(actual, expected, tolerance = 1e-10) {
  stopifnot(isTRUE(all.equal(actual, expected, tolerance = tolerance)))
}
fails <- function(expr, text) {
  error <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(error, "error"), grepl(text, conditionMessage(error), ignore.case = TRUE))
}
compare <- function(object, dt, type = NULL, mark = NULL, initial = .1, ...) {
  snapshot <- serialize(list(object, dt, type, mark, initial), NULL)
  for (infer in c(FALSE, TRUE)) {
    reference <- evaluate("r", object, dt, type, mark, initial, infer, ...)
    close(evaluate("cpp", object, dt, type, mark, initial, infer, ...), reference)
    close(evaluate("auto", object, dt, type, mark, initial, infer, ...), reference)
  }
  stopifnot(identical(snapshot, serialize(list(object, dt, type, mark, initial), NULL)))
}

# Different dimensions, unequal rates, non-square kernels, overlapping and
# reordered mappings, marks, and fixed signed coefficients.
set.seed(103)
for (dimensions in c(1L, 2L, 3L, 4L)) {
  for (kernels in unique(c(dimensions, 2L * dimensions, 5L))) {
    alpha <- matrix(runif(dimensions * kernels, .01, .08), dimensions)
    beta <- matrix(runif(dimensions * kernels, .5, 3), dimensions)
    mapping <- lapply(seq_len(dimensions), function(i) {
      unique(c(((i - 1L) %% kernels) + 1L, ((i + 1L) %% kernels) + 1L))
    })
    object <- new("hspec", mu = rep(.4, dimensions), alpha = alpha, beta = beta,
                  eta = matrix(runif(dimensions * kernels, -.01, .02), dimensions),
                  type_col_map = mapping)
    dt <- c(0, 0, rexp(38))
    type <- c(0, sample.int(dimensions, 39, replace = TRUE))
    mark <- c(0, sample(1:3, 39, replace = TRUE))
    initial <- matrix(runif(dimensions * kernels, -.03, .1), dimensions)
    compare(object, dt, type, mark, initial, N0 = 2, Nc0 = 3)
    object@beta[] <- 1.2
    compare(object, dt, type, mark, initial)
  }
}

object <- new("hspec", mu = c(.2, .3), alpha = matrix(.1, 2, 2), beta = matrix(1, 2, 2))
for (types in list(0, c(0, 1), c(0, 1, 1), c(0, 1, 2))) {
  compare(object, c(0, rep(.2, length(types) - 1L)), types)
}
compare(object, c(0, .2, .3, 0), c(0, 1, 2, 1), initial = NULL)
compare(object, c(0, .2, .3), c(1e300, 1, 2)) # Initial row is not an event.

# Extreme scales, zero intensity, negative event intensity, and overflow keep
# the R implementation's special values rather than being clipped to zero.
for (beta in c(1e-300, 1e-18, 1, 1e300)) {
  h <- new("hspec", mu = .2, alpha = .1, beta = beta)
  compare(h, c(0, 0, 1e-300, 1, 1e300))
}
for (initial in c(0, -.1, 1e308)) {
  h <- new("hspec", mu = 0, alpha = 1e308, beta = 1)
  compare(h, c(0, 0, 0, 1000), initial = initial)
}

# Baseline replay preserves R call counts, index values and index types.
called <- integer()
state <- new("hspec", mu = function(n) {
  called <<- c(called, n)
  matrix(.2 + .01 * n)
}, alpha = .1, beta = 1)
called <- integer()
reference <- evaluate("r", state, c(0, .2, .3), infer = TRUE)
reference_calls <- called
called <- integer()
close(evaluate("auto", state, c(0, .2, .3), infer = TRUE), reference)
stopifnot(identical(called, reference_calls))
called <- integer()
short_reference <- evaluate("r", state, c(0, .2))
short_calls <- called
called <- integer()
close(evaluate("cpp", state, c(0, .2)), short_reference)
stopifnot(identical(called, short_calls))
parameter_function <- new("hspec", mu = .2,
  alpha = function(param = c(a = .1)) matrix(param[1]), beta = 1)
close(evaluate("auto", parameter_function, c(0, .2)),
      evaluate("r", parameter_function, c(0, .2)))
close(evaluate("cpp", parameter_function, c(0, .2)),
      evaluate("r", parameter_function, c(0, .2)))
plain <- new("hspec", mu = .2, alpha = .1, beta = 1)
for (slot_name in c("impact", "dmark")) {
  custom <- plain
  methods::slot(custom, slot_name) <- if (slot_name == "impact") function(n) matrix(.01) else function(n) .5
  close(evaluate("auto", custom, c(0, .2, .3), infer = TRUE),
        evaluate("r", custom, c(0, .2, .3), infer = TRUE))
  if (slot_name == "impact") close(evaluate("cpp", custom, c(0, .2)),
                                   evaluate("r", custom, c(0, .2))) else
    fails(evaluate("cpp", custom, c(0, .2)), "does not support")
}
frpp <- set_residual(plain, c(rate = 2), dresidual = dexp,
                     presidual = pexp, qresidual = qexp)
close(evaluate("auto", frpp, c(0, .2, .3), infer = TRUE),
      evaluate("r", frpp, c(0, .2, .3), infer = TRUE))
fails(evaluate("cpp", frpp, c(0, .2)), "a residual distribution")
plain@rmark <- function(...) stop("Must not generate marks in logLik")
plain@rresidual <- function(...) stop("Must not generate residuals in logLik")
compare(plain, c(0, .2, .3))

# R numeric subclasses may define their own arithmetic. Auto must preserve
# those semantics instead of silently reading the underlying doubles.
check_classed_marks <- function() {
  method_names <- c("[.emhawkes_test_mark", "Ops.emhawkes_test_mark")
  assign(method_names[1], function(x, i, ...) {
    structure(NextMethod("["), class = "emhawkes_test_mark")
  }, envir = .GlobalEnv)
  assign(method_names[2], function(e1, e2) {
    value <- NextMethod()
    if (.Generic == "-") value + 1 else value
  }, envir = .GlobalEnv)
  on.exit(rm(list = method_names, envir = .GlobalEnv))
  object <- new("hspec", mu = .2, alpha = .1, beta = 1, eta = .1)
  marks <- structure(c(0, 2, 3), class = "emhawkes_test_mark")
  reference <- evaluate("r", object, c(0, .2, .3), mark = marks, infer = TRUE)
  close(evaluate("auto", object, c(0, .2, .3), mark = marks, infer = TRUE), reference)
  fails(evaluate("cpp", object, c(0, .2, .3), mark = marks), "classed observation vectors")
  stopifnot(!isTRUE(all.equal(reference,
    evaluate("cpp", object, c(0, .2, .3), mark = unclass(marks), infer = TRUE))))
}
check_classed_marks()
for (slot_name in c("mu", "type_col_map")) {
  classed <- new("hspec", mu = .2, alpha = .1, beta = 1)
  methods::slot(classed, slot_name) <- if (slot_name == "mu") {
    structure(matrix(.2), class = c("matrix", "array"))
  } else structure(list(1L), class = "list")
  close(evaluate("auto", classed, c(0, .2, .3)), evaluate("r", classed, c(0, .2, .3)))
  fails(evaluate("cpp", classed, c(0, .2)), "classed")
}
classed@type_col_map <- list(structure(1L, class = "integer"))
close(evaluate("auto", classed, c(0, .2)), evaluate("r", classed, c(0, .2)))
fails(evaluate("cpp", classed, c(0, .2)), "classed kernel mapping")

# A deterministic native likelihood neither creates nor advances R's seed.
set.seed(104)
seed <- .Random.seed
evaluate("cpp", plain, c(0, .2, .3), infer = TRUE)
stopifnot(identical(seed, .Random.seed))
rm(.Random.seed, envir = .GlobalEnv)
evaluate("cpp", plain, c(0, .2, .3))
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))

# Validate the native boundary too, so malformed internal calls cannot index
# outside allocated buffers. No blanket exception-to-minus-infinity fallback.
native <- getFromNamespace("hawkes_loglik_cpp", "emhawkes")
args <- list(dt = c(0, .2), type = c(0L, 1L), mark = c(0, 1), mu = .2,
             alpha = matrix(.1), beta = matrix(1), eta = NULL,
             initial = matrix(.1), mapping = list(0L), infer = FALSE)
for (invalid in list(list(type = c(0L, NA_integer_)), list(dt = c(1, .2)),
                     list(mu = c(.2, .3)), list(beta = matrix(0)),
                     list(eta = matrix(.1, 1, 2)), list(mapping = list(1L)),
                     list(mapping = list(c(0L, 0L))), list(mapping = list(integer())))) {
  bad_args <- args
  bad_args[names(invalid)] <- invalid
  fails(do.call(native, bad_args), "native")
}
