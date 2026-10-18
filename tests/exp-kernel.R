library(emhawkes)

fails <- function(expr, text) {
  error <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(error, "error"), grepl(text, conditionMessage(error), fixed = TRUE))
}
close <- function(a, b) stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b),
                                                   tolerance = 1e-12)))
old <- options(emhawkes.engine = "auto")
set.seed(871)
comparisons <- 0L
for (d in 1:3) for (multiple in 1:2) for (marked in c(FALSE, TRUE)) {
  k <- d * multiple
  mapping <- lapply(seq_len(d), function(i) rev(unique(c(i, k))))
  h <- new("hspec", mu = rep(.4, d), alpha = matrix(seq_len(d * k) / 500, d),
    beta = matrix(rep(seq_len(d) / 2 + 1, k), d),
    eta = if (marked) matrix(seq_len(d * k) / 900, d) else NULL,
    type_col_map = mapping)
  kernel <- as_exp_kernel(h)
  frozen <- serialize(kernel, NULL)
  for (size in c(1L, 2L, 64L)) {
    dt <- c(0, stats::rexp(size - 1L))
    if (size > 2L) dt[3L] <- 0
    type <- c(99, rep(seq_len(d), length.out = size - 1L))
    mark <- c(-99, rep(c(.5, 1, 2), length.out = size - 1L))
    initial <- matrix(seq_len(d * k) / 100, d)
    # Include genuinely negative components without changing the kernel rules.
    initial[1, 1] <- -.03
    seed <- .Random.seed
    result <- exp_kernel_replay(kernel, dt, type, mark, initial, history = TRUE)
    lean <- exp_kernel_replay(kernel, dt, type, mark, initial)
    stopifnot(identical(seed, .Random.seed), identical(frozen, serialize(kernel, NULL)),
      identical(lean, result[names(lean)]), !"component_after" %in% names(lean),
      identical(result$kernel_integral[1, ], rep(0, d)),
      identical(as.numeric(result$component_before[1, ]), as.vector(t(initial))))
    for (engine in c("r", "cpp")) {
      options(emhawkes.engine = engine)
      reference <- logLik(h, dt, type, mark, lambda_component0 = initial, infer = TRUE)
      stopifnot(identical(unname(reference$lambda_component), result$component_before),
        identical(unname(reference$rambda_component), result$component_after),
        identical(unname(reference$lambda), result$kernel_before + .4))
      comparisons <- comparisons + 1L
    }
    final <- matrix(result$component_after[size, ], d, byrow = TRUE)
    stopifnot(identical(final, result$final_component))
    # Independent all-past-events oracle, not a recurrence through the helpers.
    for (n in seq_len(size - 1L) + 1L) {
      # Integrate in local interval time. Subtracting large cumulative arrival
      # times can lose digits for short intervals and corrupt the oracle itself.
      at <- function(offset) {
        value <- initial * exp(-h@beta * (sum(dt[seq_len(n - 1L)]) + offset))
        if (n > 2L) for (j in seq.int(2L, n - 1L)) {
          jump <- matrix(0, d, k)
          columns <- mapping[[type[j]]]
          jump[, columns] <- h@alpha[, columns]
          if (marked) jump[, columns] <- jump[, columns] + h@eta[, columns] * (mark[j] - 1)
          age <- if (j == n - 1L) 0 else sum(dt[seq.int(j + 1L, n - 1L)])
          value <- value + jump * exp(-h@beta * (age + offset))
        }
        rowSums(value)
      }
      close(at(dt[n]), result$kernel_before[n, ])
      for (i in seq_len(d)) {
        integral <- stats::integrate(function(times)
          vapply(times, function(t) at(t)[i], numeric(1)),
          0, dt[n], rel.tol = 1e-12)$value
        close(integral, result$kernel_integral[n, i])
      }
    }
  }
}
cat(comparisons, "exact R/C++ state comparisons and independent integral checks passed.\n")

# Small beta * dt, exact-zero intervals and zero components must stay finite.
tiny <- exp_kernel_spec(0, .Machine$double.xmin)
z <- exp_kernel_replay(tiny, c(0, .Machine$double.xmin), lambda_component0 = 1)
stopifnot(z$kernel_integral[2, 1] == .Machine$double.xmin,
          z$final_component[1, 1] == 1)
zero <- exp_kernel_replay(tiny, c(0, 1e308), lambda_component0 = 0)
stopifnot(all(is.finite(unlist(zero))), all(unlist(zero) == 0))

# Conversion reads numeric slots without running unrelated or unsupported code.
h <- new("hspec", mu = .4, alpha = .1, beta = 1)
h@mu <- function(...) stop("baseline executed")
h@rmark <- function(...) stop("generator executed")
stopifnot(inherits(as_exp_kernel(h), "exp_kernel_spec"))
h@alpha <- function(...) stop("alpha executed")
seed <- .Random.seed
fails(as_exp_kernel(h), "alpha must")
h@alpha <- matrix(.1)
h@impact <- function(...) stop("impact executed")
fails(as_exp_kernel(h), "impact callbacks")
stopifnot(identical(seed, .Random.seed))

fails(exp_kernel_spec(c(.1, .2), 1), "matrix or scalar")
fails(exp_kernel_spec(matrix(0, 0, 1), matrix(1, 0, 1)), "matrix or scalar")
fails(exp_kernel_spec(.1, 0), "positive")
fails(exp_kernel_spec(matrix(.1, 1, 2), matrix(c(1, 2), 1), type_col_map = list(1:2)), "equal beta")
fails(exp_kernel_spec(matrix(.1, 2, 3), matrix(1, 2, 3)), "type_col_map")
fails(exp_kernel_spec(.1, 1, eta = matrix(.1, 2, 1)), "eta must match")
fails(exp_kernel_spec(.1, 1, type_col_map = list(c(1, 1))), "distinct integer")
fails(exp_kernel_spec(.1, 1, type_col_map = list(2)), "distinct integer")
k <- exp_kernel_spec(.1, 1)
fails(exp_kernel_replay(k, c(0, .1)), "explicit finite")
fails(exp_kernel_replay(k, c(0, .1), lambda_component0 = "auto"), "explicit finite")
fails(exp_kernel_replay(k, c(0, .1), lambda_component0 = NA_real_), "finite numeric")
fails(exp_kernel_replay(k, c(0, -.1), lambda_component0 = 0), "nonnegative")
fails(exp_kernel_replay(k, c(0, .1), c(0, 2), lambda_component0 = 0), "valid integer")
fails(exp_kernel_replay(k, c(0, .1), mark = 1, lambda_component0 = 0), "mark must")
fails(exp_kernel_replay(k, 0, lambda_component0 = 0, history = NA), "history must")
# Derived metadata cannot override validated matrix dimensions.
k$dimens <- 999
stopifnot(as_exp_kernel(k)$dimens == 1L)
rm(".Random.seed", envir = .GlobalEnv)
invisible(exp_kernel_replay(exp_kernel_spec(.1, 1), 0, lambda_component0 = 0))
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
options(old)
