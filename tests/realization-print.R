library(emhawkes)

h <- new("hspec", mu = .5, alpha = .1, beta = 1)
set.seed(17)
path <- hsim(h, 3, lambda_component0 = 0)
rng <- .Random.seed
old_options <- options(digits = 12)
path_display <- function(lines) {
  lines[seq.int(match("Realized path :", lines) + 1L, length(lines))]
}
for (method in list(print, summary)) {
  displayed <- path_display(capture.output(result <- withVisible(method(path, n = 1))))
  stopifnot(identical(getOption("digits"), 12L),
    identical(result$value, path), !result$visible,
    any(grepl("[1,]", displayed, fixed = TRUE)))
  displayed <- path_display(capture.output(method(path, n = 0)))
  stopifnot(!any(grepl("[1,]", displayed, fixed = TRUE)),
    any(grepl("3 more rows", displayed, fixed = TRUE)))
  displayed <- path_display(capture.output(method(path, n = Inf)))
  stopifnot(any(grepl("[3,]", displayed, fixed = TRUE)),
    !any(grepl("more rows", displayed, fixed = TRUE)))
  for (n in list(NA_real_, NaN, -1, -Inf, 1.5, numeric(), c(1, 2), "1", 1+1i)) {
    error <- tryCatch(method(path, n = n), error = identity)
    stopifnot(inherits(error, "error"), identical(getOption("digits"), 12L))
  }
  # An error after output starts must also restore the caller's options.
  broken <- path
  broken$lambda <- NULL
  displayed <- capture.output(error <- tryCatch(method(broken), error = identity))
  stopifnot(inherits(error, "error"), identical(getOption("digits"), 12L))
}
stopifnot(identical(rng, .Random.seed))
options(old_options)
