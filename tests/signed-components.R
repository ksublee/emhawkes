library(emhawkes)

close <- function(x, y) stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = 1e-9)))
fails <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))

# With equal decay rates, different component decompositions give the same path.
h <- new("hspec", mu = c(.2, .3), alpha = matrix(.1, 2, 2), beta = matrix(1, 2, 2))
signed <- matrix(c(-.5, .75, -.25, .5), 2, byrow = TRUE)
nonnegative <- matrix(c(.25, 0, .25, 0), 2, byrow = TRUE)
set.seed(31)
path <- hsim(h, size = 60, lambda_component0 = signed)
set.seed(31)
reference <- hsim(h, size = 60, lambda_component0 = nonnegative)
close(path$inter_arrival, reference$inter_arrival)
close(path$type, reference$type)
close(path$lambda, reference$lambda)
inferred <- logLik(h, path$inter_arrival, path$type, lambda_component0 = signed, infer = TRUE)
close(inferred$lambda, path$lambda)
close(logLik(h, path$inter_arrival, path$type, lambda_component0 = signed),
      logLik(h, path$inter_arrival, path$type, lambda_component0 = nonnegative))

# A state-dependent reset creates negative components while preserving positive row sums.
reset <- function(param = c(target = .1), lambda_component_n) {
  change <- param[1] - rowSums(lambda_component_n)
  matrix(rep(change, ncol(lambda_component_n)), nrow = 2)
}
h <- new("hspec", mu = c(.2, .3), alpha = matrix(0, 2, 2), beta = matrix(1, 2, 2), impact = reset)
aggregate <- new("hspec", mu = c(.2, .3), alpha = matrix(0, 2, 1), beta = matrix(1, 2, 1),
                 type_col_map = list(1, 1), impact = reset)
set.seed(32)
path <- hsim(h, size = 60, lambda_component0 = .5)
set.seed(32)
reference <- hsim(aggregate, size = 60, lambda_component0 = 1)
stopifnot(any(path$rambda_component < 0))
close(path$inter_arrival, reference$inter_arrival)
close(path$type, reference$type)
close(path$lambda, reference$lambda)
inferred <- logLik(h, path$inter_arrival, path$type, lambda_component0 = .5, infer = TRUE)
close(inferred$lambda, path$lambda)
close(inferred$rambda, path$rambda)

# Signed fast kernels in the synchronized model preserve decreasing excitation.
alpha <- matrix(c(4, -1, 4, 1, 4, 1, 4, -1), 2, byrow = TRUE)
beta <- matrix(rep(c(10.56, 13.44), 4), 2, byrow = TRUE)
initial <- matrix(1, 2, 4)
synchronized <- new("hspec", mu = c(.2, .3), alpha = alpha, beta = beta,
                    type_col_map = list(1:2, 3:4))
set.seed(33)
path <- hsim(synchronized, size = 60, lambda_component0 = initial)
stopifnot(any(path$rambda_component < 0))

# Evaluate intensity from the full event history and numerically integrate it,
# independently of the package's recursive intensity and integral calculations.
times <- cumsum(path$inter_arrival)
intensity <- function(t, i, n) {
  value <- synchronized@mu[i] + sum(initial[i, ] * exp(-beta[i, ] * t))
  if (n > 2L) for (j in 2:(n - 1L)) {
    columns <- synchronized@type_col_map[[path$type[j]]]
    value <- value + sum(alpha[i, columns] * exp(-beta[i, columns] * (t - times[j])))
  }
  value
}
reference <- path$lambda
integral <- 0
for (n in 2:length(times)) for (i in 1:2) {
  reference[n, i] <- intensity(times[n], i, n)
  integral <- integral + integrate(function(t) vapply(t, intensity, numeric(1), i = i, n = n),
                                   times[n - 1L], times[n], rel.tol = 1e-10)$value
}
close(path$lambda, reference)
expected <- sum(log(reference[cbind(2:length(times), path$type[-1])])) - integral
close(logLik(synchronized, path$inter_arrival, path$type, lambda_component0 = initial), expected)
inferred <- logLik(synchronized, path$inter_arrival, path$type, lambda_component0 = initial, infer = TRUE)
close(inferred$lambda, reference)
close(inferred$loglikelihood, expected)

# Rate groups, rather than the column order, determine admissibility.
permutation <- c(2, 1, 4, 3)
permuted <- new("hspec", mu = synchronized@mu, alpha = alpha[, permutation],
               beta = beta[, permutation], type_col_map = list(1:2, 3:4))
set.seed(33)
reference <- hsim(permuted, size = 60, lambda_component0 = initial[, permutation])
close(path$inter_arrival, reference$inter_arrival)
close(path$type, reference$type)
close(path$lambda, reference$lambda)

# Negative excitation and signed kernels that fail the monotonicity condition are rejected.
fails(hsim(h, size = 2, lambda_component0 = matrix(-.1, 2, 2)))
unequal <- new("hspec", mu = c(.2, .3), alpha = matrix(.1, 2, 2), beta = matrix(c(1, 2, 1, 2), 2, byrow = TRUE))
fails(hsim(unequal, size = 2, lambda_component0 = signed))
