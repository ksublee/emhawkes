# Before changing the package, install the R baseline into a separate library.
# R_LIBS=<baseline-lib> Rscript tools/check-rcpp.R snapshot <baseline.rds>
# R_LIBS=<native-lib>   Rscript tools/check-rcpp.R compare <baseline.rds> <report.csv>
library(emhawkes)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || !args[1] %in% c("snapshot", "compare")) {
  stop("Use: check-rcpp.R snapshot|compare baseline.rds [report.csv]")
}

collect <- function(engine = "r") {
  previous <- options(emhawkes.engine = engine)
  on.exit(options(previous))
  set.seed(731)
  dt <- c(0, 0, rexp(298, rate = 2))
  types <- c(0, sample(1:2, 299, replace = TRUE))
  marks <- c(0, sample(1:3, 299, replace = TRUE))
  alpha <- matrix(c(.12, .04, .06, .15), 2)
  beta <- matrix(c(.8, 1.2, 1.6, 2), 2)
  ordinary <- new("hspec", mu = c(.2, .3), alpha = alpha, beta = beta)
  marked <- ordinary
  marked@eta <- matrix(.02, 2, 2)
  multi <- new("hspec", mu = c(.2, .3), alpha = cbind(alpha, alpha / 2),
               beta = cbind(beta, beta * 2), type_col_map = list(c(3, 1), c(4, 2)))
  signed <- ordinary
  signed@alpha[1, 2] <- -.015
  specifications <- list(ordinary = ordinary, marked = marked, multi = multi, signed = signed)
  likelihoods <- lapply(specifications, function(object) {
    list(scalar = logLik(object, dt, types, marks, lambda_component0 = .1),
         inference = logLik(object, dt, types, marks, lambda_component0 = .1, infer = TRUE))
  })
  likelihoods$default_initial <- suppressMessages(logLik(ordinary, dt, types, marks, infer = TRUE))

  # Real simulated data, with numeric parameters estimated from the same start.
  set.seed(732)
  truth <- new("hspec", mu = .35, alpha = .3, beta = 1.5)
  sample <- hsim(truth, size = 1001, lambda_component0 = .1)
  start <- new("hspec", mu = .25, alpha = .2, beta = 1.2)
  fit <- hfit(start, sample$inter_arrival, lambda_component0 = .1,
    constraints = list(ineqA = diag(3), ineqB = rep(-1e-6, 3)),
    control = list(iterlim = 200, reltol = 1e-9))
  set.seed(733)
  truth <- new("hspec", mu = c(.25, .35), alpha = matrix(.06, 2, 4),
               beta = matrix(rep(c(.8, 1.2, 1.6, 2), each = 2), 2),
               type_col_map = list(c(1, 3), c(2, 4)))
  sample <- hsim(truth, size = 1001, lambda_component0 = .1)
  start <- truth
  start@mu <- matrix(c(.2, .3), 2)
  start@alpha[] <- .04
  # Keep four distinct decay parameters fixed; estimate two baselines and a
  # shared excitation parameter. Multi-kernel sharing is still exercised.
  multi_fit <- suppressMessages(hfit(start, sample$inter_arrival, sample$type,
    lambda_component0 = .1, fixed = c(rep(FALSE, 3), rep(TRUE, 4)),
    constraints = list(ineqA = diag(7), ineqB = rep(-1e-6, 7)),
    control = list(iterlim = 200, reltol = 1e-9)))
  fit_summary <- function(object) list(coef = coef(object), logLik = as.numeric(logLik(object)),
    code = object$code, gradient = object$gradient, hessian = object$hessian)
  list(likelihoods = likelihoods, fits = list(univariate = fit_summary(fit), multi = fit_summary(multi_fit)))
}

if (args[1] == "snapshot") {
  saveRDS(list(results = collect(), session = sessionInfo()), args[2])
  cat("Saved pre-Rcpp results to", args[2], "\n")
} else {
  baseline <- readRDS(args[2])$results
  reports <- list()
  for (engine in c("r", "cpp", "auto")) {
    result <- collect(engine)
    stopifnot(isTRUE(all.equal(result$likelihoods, baseline$likelihoods, tolerance = 1e-10)))
    max_difference <- function(x, y) max(abs(unlist(x) - unlist(y)))
    for (name in names(result$fits)) {
      current <- result$fits[[name]]
      previous <- baseline$fits[[name]]
      stopifnot(abs(current$logLik - previous$logLik) < 1e-7,
                identical(names(current$coef), names(previous$coef)),
                max_difference(current$coef, previous$coef) < 1e-4,
                current$code == 0, previous$code == 0)
      reports[[length(reports) + 1L]] <- data.frame(engine = engine, fit = name,
        likelihood_max_abs = max_difference(result$likelihoods, baseline$likelihoods),
        fit_loglik_abs = abs(current$logLik - previous$logLik),
        coefficient_max_abs = max_difference(current$coef, previous$coef),
        gradient_max_abs = max_difference(current$gradient, previous$gradient),
        hessian_max_abs = max_difference(current$hessian, previous$hessian),
        convergence_code = current$code)
    }
  }
  reports <- do.call(rbind, reports)
  print(reports, row.names = FALSE)
  if (length(args) >= 3L) write.csv(reports, args[3], row.names = FALSE)
}
