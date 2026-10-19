library(emhawkes)
fails <- function(expr, text) {
  error <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(error, "error"), grepl(text, conditionMessage(error), fixed = TRUE))
}
set.seed(984)
saved <- .Random.seed
comparisons <- 0L
for (d in 1:3) for (multiple in 1:2) for (marked in c(FALSE, TRUE)) {
  k <- d * multiple
  kernel <- exp_kernel_spec(matrix(seq_len(d * k) / 100, d),
    matrix(rep(seq_len(d) / 2 + 1, k), d),
    eta = if (marked) matrix(.02, d, k) else NULL,
    type_col_map = lapply(seq_len(d), function(i) rev(unique(c(i, k)))))
  component <- initial <- matrix(seq_len(d * k) / 20, d)
  initial[1, 1] <- component[1, 1] <- -.02
  dt <- c(0, .2, 0, .5)
  types <- c(99, rep(seq_len(d), length.out = 3))
  marks <- c(0, .5, 2, 1.5)
  reference <- exp_kernel_replay(kernel, dt, types, marks, initial, history = TRUE)
  frozen <- serialize(kernel, NULL)
  for (n in 2:length(dt)) {
    step <- exp_kernel_step(kernel, component, dt[n], types[n], marks[n])
    lean <- exp_kernel_step(kernel, component, dt[n], types[n], marks[n], integrate = FALSE)
    stopifnot(identical(step[names(lean)], lean),
      identical(as.vector(t(step$before)), reference$component_before[n, ]),
      identical(as.vector(t(step$after)), reference$component_after[n, ]),
      identical(rowSums(step$integral), reference$kernel_integral[n, ]))
    component <- step$after
    comparisons <- comparisons + 1L
  }
  stopifnot(identical(component, reference$final_component),
    identical(frozen, serialize(kernel, NULL)), identical(saved, .Random.seed))
}
kernel <- exp_kernel_spec(matrix(.1, 2, 2), matrix(c(1, 2, 1, 2), 2))
mu <- c(.2, .4)
stopifnot(identical(exp_kernel_wait(kernel, 0, mu, c(.3, .7)), c(.3, .7) / mu),
  identical(exp_kernel_wait(kernel, 0, 0, c(0, 1)), c(0, Inf)),
  identical(exp_kernel_wait(kernel, .1, 0, c(Inf, .1)), c(Inf, Inf)))
component <- matrix(c(-.02, .1, .12, .2), 2)
target <- c(.4, 1.3)
wait <- exp_kernel_wait(kernel, component, mu, target)
for (i in 1:2) {
  hazard <- integrate(function(t) mu[i] + sum(component[i, ]) *
    exp(-kernel$beta[i, 1] * t), 0, wait[i], rel.tol = 1e-12)$value
  stopifnot(isTRUE(all.equal(hazard, target[i], tolerance = 1e-11)))
}
stopifnot(identical(exp_kernel_step(kernel, component, 0, integrate = FALSE)$after, component),
  identical(saved, .Random.seed))
fails(exp_kernel_step(kernel, 0, -1), "dt")
fails(exp_kernel_step(kernel, 0, 0, 3), "type")
fails(exp_kernel_step(kernel, 0, 0, 1, NA_real_), "mark")
fails(exp_kernel_step(kernel, numeric(0), 0), "component")
fails(exp_kernel_step(kernel, matrix(0, 1, 4), 0), "component")
fails(exp_kernel_step(kernel, 0, 0, integrate = NA), "integrate")
fails(exp_kernel_wait(kernel, 0, -1, c(1, 1)), "mu")
fails(exp_kernel_wait(kernel, 0, 1, 1), "residual")
fails(exp_kernel_wait(kernel, 0, 1, c(NA, 1)), "residual")
fails(exp_kernel_wait(kernel, -1, 1, c(1, 1)), "nonnegative")
stopifnot(identical(saved, .Random.seed))
cat(comparisons, "event-wise state comparisons; deterministic inversion, quadrature and failure contracts passed.\n")
