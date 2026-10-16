library(emhawkes)
r_sim <- getFromNamespace("simulate_hawkes_r", "emhawkes")
data_only <- function(x) x[setdiff(names(x), "hspec")]
same <- function(h, size, initial, seed, N0 = c(2), Nc0 = c(5)) {
  # Function environments can force existing promises during evaluation.
  numeric_model <- !any(vapply(methods::slotNames(h), function(name)
    is.function(methods::slot(h,name)),logical(1)))
  before_object <- if(numeric_model) serialize(h, NULL) else NULL
  capture <- function(expr) {
    warnings <- character()
    value <- withCallingHandlers(force(expr),warning=function(w) {
      warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
    })
    list(value=value,warnings=warnings)
  }
  set.seed(seed)
  reference <- capture(r_sim(h, size, lambda_component0 = initial, N0 = N0, Nc0 = Nc0))
  rng <- .Random.seed
  set.seed(seed)
  native <- capture(hsim(h, size, lambda_component0 = initial, N0 = N0, Nc0 = Nc0))
  stopifnot(identical(data_only(reference$value), data_only(native$value)),
    identical(rng, .Random.seed),identical(reference$warnings,native$warnings))
  if(numeric_model) stopifnot(identical(before_object, serialize(h, NULL)))
  native$value
}
models <- list()
for (d in c(1L, 2L, 3L, 4L, 8L)) for (k in unique(c(1L, d, 2L*d, 5L))) {
  mapping <- if (k == d) NULL else lapply(seq_len(d), function(i) rev(unique(c(i %% k + 1L, k))))
  h <- new("hspec", mu = seq(.3, .7, length.out = d), alpha = matrix(.06 / k, d, k),
    beta = matrix(rep(seq(.8, 1.9, length.out = d), k), d, k), type_col_map = mapping)
  models[[paste(d, k, sep = ":")]] <- h
}
models$signed <- new("hspec", mu = c(.2, .3),
  alpha = matrix(c(.08,-.04,.05,-.02,.05,-.02,.08,-.04),2,byrow=TRUE),
  beta = matrix(c(1.1,1.5),2,4), type_col_map = list(c(2L,1L),c(4L,3L)))
models$poisson <- new("hspec", mu=.5,alpha=0,beta=1)
models$slow <- new("hspec",mu=1e-6,alpha=2e-7,beta=1e-5)
models$fast <- new("hspec",mu=1e5,alpha=1e3,beta=1e6)
saved_rng <- RNGkind()
for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG")) {
  RNGkind(kind, "Inversion", "Rejection")
  for (h in models) for (size in c(1L,2L,31L,251L)) for (seed in c(1L,62L,7361L)) {
    initial <- if (identical(h,models$signed)) h@alpha/2 else .03/ncol(h@beta)
    same(h,size,initial,seed,N0=seq_len(h@dimens),Nc0=seq_len(h@dimens)+2)
  }
}
do.call(RNGkind, as.list(saved_rng))

# Default initialization and legacy scalar/matrix normalization remain in R.
h <- models[["2:2"]]
for (initial in list(NULL, .04, matrix(1:4/100,2))) same(h,101,initial,511)
# Uniforms are still consumed when the excitation sum is exactly zero.
invisible(same(models$poisson,101,0,513))
# Reciprocal overflow and subnormal waits must retain R's trimming/warnings.
for (level in c(1e-320,1e-300,1e308))
  same(new("hspec",mu=level,alpha=0,beta=1),3,0,514)
# Likelihood engine options have no role in simulation dispatch.
old_options <- options(emhawkes.engine = "r")
for (engine in c("r","auto","cpp","unrelated")) {
  options(emhawkes.engine=engine)
  same(h,31,.03,512)
}
options(old_options)

# A one-row native request must neither create nor advance the RNG state.
if (exists(".Random.seed",.GlobalEnv)) rm(".Random.seed",envir=.GlobalEnv)
one <- hsim(h,1,lambda_component0=.03)
stopifnot(!exists(".Random.seed",.GlobalEnv))
set.seed(1); rng <- .Random.seed
stopifnot(identical(one,hsim(h,1,lambda_component0=.03)),identical(rng,.Random.seed))

# Reconstruct generated pre/post intensities through two independent engines.
for (h in models) {
  initial <- if (identical(h,models$signed)) h@alpha/2 else .03/ncol(h@beta)
  set.seed(913); path <- hsim(h,301,lambda_component0=initial,N0=2,Nc0=5);rng <- .Random.seed
  for (engine in c("r","cpp")) {
    op <- options(emhawkes.engine=engine)
    replay <- logLik(h,path$inter_arrival,path$type,path$mark,lambda_component0=initial,
      N0=2,Nc0=5,infer=TRUE)
    options(op)
    for (field in c("lambda","rambda","lambda_component","rambda_component"))
      stopifnot(isTRUE(all.equal(as.numeric(path[[field]]),as.numeric(replay[[field]]),tolerance=1e-12)))
    stopifnot(identical(rng,.Random.seed))
  }
  for (i in seq_len(h@dimens)) stopifnot(
    identical(path$N[,i],2+cumsum(path$type==i)),
    identical(path$Nc[,i],5+cumsum((path$type==i)*path$mark)))
}

# All unsupported specifications are classified before the native RNG scope.
unsupported <- list(
  new("hspec",mu=0,alpha=.1,beta=1),
  new("hspec",mu=c(0,.4),alpha=matrix(.01,2,2),beta=matrix(1,2,2)),
  new("hspec",mu=.4,alpha=matrix(.01,1,2),beta=matrix(c(1,2),1),type_col_map=list(1:2)),
  new("hspec",mu=.4,alpha=.1,beta=1,eta=.01),
  new("hspec",mu=function(param=c(level=.4))matrix(param[1]),alpha=.1,beta=1),
  new("hspec",mu=.4,alpha=.1,beta=1,impact=function(lambda_component_n) .01*lambda_component_n),
  new("hspec",mu=.4,alpha=.1,beta=1,rmark=function()rgeom(1,.7)+1))
unsupported[[length(unsupported)+1L]] <- set_residual(new("hspec",mu=.4,alpha=.1,beta=1),
  param=c(shape=2),rresidual=function(n,shape)rgamma(n,shape),
  dresidual=function(x,shape)dgamma(x,shape),presidual=function(q,shape)pgamma(q,shape),
  qresidual=function(p,shape)qgamma(p,shape))
trace("hawkes_simulate_cpp", where=asNamespace("emhawkes"),
  tracer=quote(stop("unsupported model reached native simulation")),print=FALSE)
for (h in unsupported) same(h,31,.03,715)
untrace("hawkes_simulate_cpp",where=asNamespace("emhawkes"))

# Initial validation failures do not consume RNG; callbacks remain covered by
# simulation-dispatch.R, including snapshots retained until after simulation.
set.seed(18);rng <- .Random.seed
error <- tryCatch(hsim(models[[1]],31,lambda_component0=-1),error=identity)
stopifnot(inherits(error,"error"),identical(rng,.Random.seed))

# An error after a draw must retain that consumed RNG state, without restart.
h <- new("hspec",mu=.4,alpha=-.2,beta=1)
set.seed(19);old_error <- tryCatch(r_sim(h,5,lambda_component0=0),error=identity);rng <- .Random.seed
set.seed(19);new_error <- tryCatch(hsim(h,5,lambda_component0=0),error=identity)
stopifnot(inherits(old_error,"error"),inherits(new_error,"error"),
  identical(conditionMessage(old_error),conditionMessage(new_error)),identical(rng,.Random.seed))
