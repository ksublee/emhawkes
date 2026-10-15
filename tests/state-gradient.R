library(emhawkes)
# JIT compilation changes serialized closures even when their values are intact.
previous_jit <- compiler::enableJIT(0)
normalize <- emhawkes:::normalize_observation
build <- emhawkes:::build_parameter_map
objective <- emhawkes:::make_hawkes_objective
restore <- emhawkes:::restore_parameters
close <- function(a, b, tol = 4e-7) stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol)))
fails <- function(expr, pattern) {
  e <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(e, 'error'), grepl(pattern, conditionMessage(e), ignore.case = TRUE))
}
central <- function(f, p, step) vapply(seq_along(p), function(k) {
  a <- b <- p; h <- step * max(1, abs(p[k])); a[k] <- a[k] + h; b[k] <- b[k] - h
  (f(a) - f(b)) / (2 * h)
}, 0)
old <- options(emhawkes.engine = 'auto')
dt <- c(0, .3, 0, .7, .2, 1.1, .05)
type <- c(0, 1, 2, 2, 1, 2, 1)
mark <- c(0, 1, 3, .5, 2, 1, 2.5)
# Baseline depends on observed counts and resolved alpha, sharing the name 'a'.
mu <- function(param = c(m = -.8, a = -.9), n, N, alpha, ...) {
  matrix(c(exp(param['m']) * (1 + N[max(1, n - 1), 1] / 30) + alpha[1] / 10,
    exp(param['a']) * (1 + n / 100)), 2)
}
alpha <- function(param = c(a = -.9)) matrix(c(exp(param), .04, .06, exp(param) / 2), 2)
beta <- function(param = c(b = .3)) matrix(c(1, 2, 3, 4) * exp(param), 2)
eta <- matrix(c(.01, .02, 0, .03), 2)
impact <- function(param = c(q = .02), n, Nc, alpha, lambda_component_n, ...) {
  matrix(param['q'] * (1 + sum(Nc[n, ]) / 50) + as.vector(alpha) / 30 +
    .02 * as.vector(lambda_component_n)^2 + .01 * sum(lambda_component_n), 2)
}
j <- list(
  mu = function(param, n, N, ...) {
    J <- matrix(0, 2, length(param), dimnames = list(NULL, names(param)))
    J[1, 'm'] <- exp(param['m']) * (1 + N[max(1, n - 1), 1] / 30)
    J[1, 'a'] <- exp(param['a']) / 10
    J[2, 'a'] <- exp(param['a']) * (1 + n / 100)
    # Every callback result may return a different column order.
    J[, if (n %% 2) rev(names(param)) else names(param), drop = FALSE]
  },
  alpha = function(param) matrix(c(exp(param), 0, 0, exp(param) / 2), dimnames = list(NULL, 'a')),
  beta = function(param) matrix(c(1, 2, 3, 4) * exp(param), dimnames = list(NULL, 'b')),
  impact = list(
    param = function(param, n, Nc, ...) {
      J <- matrix(0, 4, length(param), dimnames = list(NULL, names(param)))
      J[, 'q'] <- 1 + sum(Nc[n, ]) / 50
      J[, 'a'] <- c(exp(param['a']), 0, 0, exp(param['a']) / 2) / 30
      J[, rev(names(param)), drop = FALSE]
    },
    state = function(lambda_component_n) diag(.04 * as.vector(lambda_component_n)) + .01))
h <- new('hspec', mu = mu, alpha = alpha, beta = beta, eta = eta, impact = impact)
check <- function(model, jac, initial, dt0 = dt, type0 = type, mark0 = mark) {
  snapshot <- serialize(model, NULL)
  for (reduced in c(TRUE, FALSE)) {
    o <- normalize(model, dt0, type0, mark0, N0 = 2, Nc0 = 3)
    m <- build(model, reduced)
    f <- objective(model, o, initial, m, jacobian = jac)
    numeric <- objective(model, o, initial, m, derivatives = FALSE)
    reference <- function(p) {
      candidate <- model
      v <- restore(m, p)
      for (name in c('mu', 'alpha', 'beta', 'eta', 'impact')) slot(candidate, name) <- v[[name]]
      op <- options(emhawkes.engine = 'r'); on.exit(options(op))
      as.numeric(suppressMessages(logLik(candidate, dt0, type0, mark0,
        lambda_component0 = initial, N0 = 2, Nc0 = 3)))
    }
    p <- m$start
    direction <- rep(c(.1, -.2, .15), length.out = length(p))
    for (point in list(p, p + .02 * direction)) {
      close(f$loglik(point), numeric$loglik(point), 1e-13)
      close(f$loglik(point), reference(point), 1e-12)
      g <- f$gradient(point)
      stopifnot(identical(names(g), names(p)))
      for (step in c(1e-4, 2e-5)) close(g, central(reference, point, step))
      close(sum(g * direction), (reference(point + 2e-5 * direction) -
        reference(point - 2e-5 * direction)) / 4e-5)
      close(f$gradient(rev(point)), g, 0)
    }
    bad <- p; bad[1] <- Inf
    stopifnot(identical(as.numeric(f$loglik(bad)), -Inf), all(is.na(f$gradient(bad))))
  }
  stopifnot(identical(snapshot, serialize(model, NULL)))
}
check(h, j, matrix(c(.1, .02, .03, .15), 2))
check(h, j, NULL)
multi <- h
multi@type_col_map <- list(c(2, 1), 2)
check(multi, j, .1)
# A fixed observed baseline has no local parameters, but still needs a total J.
fixed <- new('hspec', mu = function(n) matrix(.4 + n / 100), alpha = .2, beta = 1.3)
m <- build(fixed, TRUE)
zero <- matrix(0, 1, length(m$start), dimnames = list(NULL, names(m$start)))
check(fixed, list(mu = zero), NULL, type0 = rep(1, length(dt)))
# One-dimensional state feedback, constant state Jacobian, automatic initializer.
u <- new('hspec', mu = .4, alpha = .2, beta = 1.3,
  impact = function(param = c(q = .02), lambda_component_n) matrix(param['q'] + .03 * lambda_component_n))
m <- build(u, TRUE)
pj <- matrix(0, 1, length(m$start), dimnames = list(NULL, names(m$start))); pj[, 'q'] <- 1
uj <- list(impact = list(param = pj, state = matrix(.03)))
check(u, uj, NULL, type0 = rep(1, length(dt)))
check(u, uj, .2, dt0 = c(0, 0, 1e-10, 100, .4, .6, .1), type0 = rep(1, length(dt)))
# Extreme rates use a scaled local beta parameter so independent differences
# stay in-domain and exercise the integral derivative without huge perturbations.
for(rate in c(1e-306,1e306)) {
  extreme <- u
  extreme@beta <- local({scale<-rate;function(param=c(b=0))matrix(scale*exp(param))})
  ej <- list(beta=local({scale<-rate;function(param)matrix(scale*exp(param),dimnames=list(NULL,'b'))}),
    impact=list(param=function(param) {
      J<-matrix(0,1,length(param),dimnames=list(NULL,names(param)));J[,'q']<-1;J
    },state=matrix(.03)))
  intervals <- if(rate>1) c(0,.3,.4,0,1,.5,.2)/rate else dt
  check(extreme,ej,.1,dt0=intervals,type0=rep(1,length(dt)))
}
# Three kernels, overlapping selected columns, cross-component state feedback.
k <- new('hspec', mu = c(.4, .5), alpha = matrix(.03, 2, 3),
  beta = matrix(c(.7, 1.1, 1.4, 2, 3, 4), 2), eta = matrix(.01, 2, 3),
  type_col_map = list(c(3, 1), c(2, 3)),
  impact = function(param = c(q = .01), lambda_component_n) matrix(param + .02 * sum(lambda_component_n), 2, 3))
m <- build(k, TRUE)
pj <- matrix(0, 6, length(m$start), dimnames = list(NULL, names(m$start))); pj[, 'q'] <- 1
check(k, list(impact = list(param = function(param, ...) {
  J <- matrix(0, 6, length(param), dimnames = list(NULL, names(param)));J[, 'q'] <- 1;J
}, state = matrix(.02, 6, 6))), .1)
# Observed-only impact in combination with an observed baseline.
o <- fixed; o@impact <- function(n, alpha) alpha * n / 100
m <- build(o, TRUE)
ij <- matrix(0, 1, length(m$start), dimnames = list(NULL, names(m$start)))
check(o, list(mu = zero, impact = function(param, n) {
  J <- matrix(0, 1, length(param), dimnames = list(NULL, names(param))); J[, 'alpha1'] <- n / 100; J
}), .1, type0 = rep(1, length(dt)))
# Independent Poisson formula (observed baseline at the next interval's index).
p <- new('hspec', mu = function(param = c(m = .6), n) matrix(param * (1 + n / 10)),
  alpha = function(param = numeric(0)) matrix(0), beta = 1)
m <- build(p, TRUE); po <- normalize(p, dt)
pj <- function(param, n) {J <- matrix(0, 1, length(param), dimnames=list(NULL,names(param)));J[, 'm'] <- 1+n/10;J}
f <- objective(p, po, 0, m, jacobian = list(mu = pj))
expected <- (length(dt)-1) / m$start['m'] - sum(dt[-1] * (1 + seq_along(dt)[-1] / 10))
close(f$gradient(m$start)['m'], expected, 1e-12)
# Empty paths neither call a terminal Jacobian nor an impact Jacobian.
empty_j <- j; empty_j$mu <- function(param, n, ...) {
  stopifnot(n == 1L); matrix(0, 2, length(param), dimnames=list(NULL,names(param)))
}; empty_j$impact$param <- function(...) stop('unused impact derivative')
f <- objective(h, normalize(h, 0, 0, 0), .1, build(h, TRUE), jacobian = empty_j)
close(f$loglik(build(h, TRUE)$start), 0, 0); close(f$gradient(build(h, TRUE)$start), rep(0,length(build(h,TRUE)$start)), 0)
# Explicit capability and shape errors; no silent numeric fallback.
make <- function(model = h, jac = j) objective(model, normalize(model, dt, type, mark), .1, build(model, TRUE), jacobian = jac)
bad <- j;bad$impact <- NULL;fails(make(jac=bad), 'requires jacobian.impact|requires jacobian\\$impact')
bad <- j;bad$mu <- NULL;fails(make(jac=bad), 'requires jacobian.mu')
bad <- j;bad$impact$state <- matrix(0, 3, 3);fails(make(jac=bad), '4 rows.*4 columns')
bad <- j;bad$impact$param <- function(param) matrix(0,4,length(param));fails(make(jac=bad)$gradient(build(h,TRUE)$start), 'column names')
bad <- j;bad$impact$state <- function(...) matrix(NA_real_,4,4);fails(make(jac=bad)$gradient(build(h,TRUE)$start), 'finite')
bad <- h;bad@impact <- function(lambda_component) matrix(0,2,2);fails(make(bad), 'full intensity histories')
bad <- j;bad$mu <- function(param,lambda) matrix(0,2,length(param));fails(make(jac=bad), 'latent intensity')
bad <- h;bad@impact <- function(...) {list(...)$lambda; matrix(0,2,2)}
bj <- j;bj$impact <- j$impact$param
fails(make(bad,bj)$loglik(build(bad,TRUE)$start), 'full intensity histories')
bad <- h;bad@impact <- function(lambda_component_n, ...) {list(...)$lambda;matrix(0,2,2)}
fails(make(bad)$loglik(build(bad,TRUE)$start), 'full intensity histories')
# Value and derivative passes must replay exactly the same deterministic model.
changing <- u
calls <- 0L
changing@impact <- function(param=c(q=.02),lambda_component_n) {
  calls <<- calls+1L
  matrix(param['q']+calls/100)
}
f <- objective(changing,normalize(changing,dt),.1,build(changing,TRUE),jacobian=uj)
fails(f$gradient(build(changing,TRUE)$start),'changed between value and gradient')
# A singular automatic initializer requires explicit initial components.
singular <- u; singular@alpha <- singular@beta
fails(objective(singular, normalize(singular, dt), NULL, build(singular,TRUE),
  jacobian=uj)$loglik(build(singular,TRUE)$start), 'singular automatic')
# The nested derivative contract is retained in a fitted current-state model.
state_fit <- suppressMessages(hfit(u, dt, lambda_component0=.1, method='NR',
  fixed=c(FALSE,TRUE,TRUE,TRUE), gradient='analytic', jacobian=uj))
state_num <- suppressMessages(hfit(u, dt, lambda_component0=.1, method='NR',
  fixed=c(FALSE,TRUE,TRUE,TRUE)))
close(coef(state_fit),coef(state_num),1e-5)
close(state_fit$maximum,state_num$maximum,1e-9)
stopifnot(identical(attr(state_fit,'emhawkes')$jacobian,uj),
  fit_diagnostics(state_fit)$computation$gradient$used=='analytic')
# A callback's current-state input is also isolated in derivative calls.
state_seen <- list()
snapshot_j <- uj
snapshot_j$impact$state <- function(lambda_component_n) {
  state_seen[[length(state_seen)+1L]] <<- local({x<-lambda_component_n;function() x})
  lambda_component_n[] <- 777
  matrix(.03)
}
f <- objective(u, normalize(u, dt), .1, build(u,TRUE), jacobian=snapshot_j)
close(f$gradient(build(u,TRUE)$start),central(
  objective(u,normalize(u,dt),.1,build(u,TRUE),derivatives=FALSE)$loglik,build(u,TRUE)$start,2e-5))
stopifnot(all(vapply(state_seen,function(x)all(x()<777),TRUE)))
# Retained snapshots survive call reuse, GC, and later callback mutation.
saved <- list()
u@impact <- function(param = c(q = .02), lambda_component_n) {
  saved[[length(saved)+1L]] <<- local({x<-lambda_component_n;function() x})
  gc(FALSE); lambda_component_n[] <- 999; matrix(param['q'])
}
uj$impact$state <- matrix(0)
f <- objective(u, normalize(u, dt), .1, build(u, TRUE), jacobian = uj)
invisible(f$loglik(build(u, TRUE)$start)); first <- lapply(saved,function(x)x())
invisible(f$gradient(build(u, TRUE)$start))
stopifnot(identical(lapply(saved[seq_along(first)],function(x)x()), first), all(vapply(first,function(x)all(x<999),TRUE)))
# Fit reconstruction and diagnostics preserve the new nested derivative contract.
a <- list(object=p,inter_arrival=c(0,.3,.7,.2,.5,.9),lambda_component0=0,method='NR',fixed=2)
fit <- suppressMessages(do.call(hfit,c(a,list(gradient='analytic',jacobian=list(mu=pj)))))
num <- suppressMessages(do.call(hfit,a))
close(coef(fit),coef(num),1e-5);close(as.numeric(logLik(fit)),as.numeric(logLik(num)),1e-9)
stopifnot(fit_diagnostics(fit)$computation$gradient$used == 'analytic', identical(attr(fit,'emhawkes')$jacobian,list(mu=pj)))
close(as.numeric(logLik(fitted_hspec(fit),a$inter_arrival,lambda_component0=0)),as.numeric(logLik(fit)),1e-10)
options(old)
invisible(compiler::enableJIT(previous_jit))
