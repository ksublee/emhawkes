library(emhawkes)
r_sim <- getFromNamespace("simulate_hawkes_r","emhawkes")
compare <- function(h,size,initial,seed) {
  before <- serialize(h,NULL)
  capture <- function(expr) {
    warnings <- character()
    value <- withCallingHandlers(tryCatch(force(expr),error=identity),warning=function(w) {
      warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
    })
    list(value=value,warnings=warnings)
  }
  set.seed(seed)
  ref <- capture(r_sim(h,size,lambda_component0=initial,N0=2,Nc0=5));rng <- .Random.seed
  set.seed(seed)
  out <- capture(hsim(h,size,lambda_component0=initial,N0=2,Nc0=5))
  stopifnot(identical(rng,.Random.seed),identical(ref$warnings,out$warnings),
    identical(before,serialize(h,NULL)))
  if(inherits(ref$value,"error")) {
    stopifnot(inherits(out$value,"error"),identical(conditionMessage(ref$value),conditionMessage(out$value)))
  } else stopifnot(identical(ref$value,out$value))
  out$value
}
models <- list()
for(d in c(1L,2L,4L,8L)) for(k in unique(c(2L,d,2L*d,5L))) {
  mapping <- if(k==d) NULL else lapply(seq_len(d),function(i) rev(unique(c(i%%k+1L,k))))
  beta <- outer(seq(.7,1.4,length.out=d),seq(.5,2,length.out=k))
  # One equal-rate row must also use thinning when another row is unequal.
  if(d>1) beta[1,] <- 1
  models[[paste(d,k,sep=":")]] <- new("hspec",mu=seq(.1,.5,length.out=d),
    alpha=matrix(.1/k,d,k),beta=beta,type_col_map=mapping)
}
models$signed <- new("hspec",mu=c(.2,.3),
  alpha=matrix(c(.08,-.02,.08,-.02),2,4,byrow=TRUE),
  beta=matrix(c(1,2,1,2),2,4,byrow=TRUE),type_col_map=list(c(2L,1L),c(4L,3L)))
permutation <- c(4L,1L,3L,2L)
models$permuted <- new("hspec",mu=models$signed@mu,
  alpha=models$signed@alpha[,permutation],beta=models$signed@beta[,permutation],
  type_col_map=list(c(4L,2L),c(1L,3L)))
models$rejection <- new("hspec",mu=c(.01,.02),alpha=matrix(.1,2,4),
  beta=matrix(seq(10,80,length.out=8),2,4),type_col_map=list(1:2,3:4))
saved <- RNGkind()
for(kind in c("Mersenne-Twister","L'Ecuyer-CMRG")) {
  RNGkind(kind,"Inversion","Rejection")
  for(h in models) for(size in c(1L,2L,31L,251L)) for(seed in c(1L,62L,7361L)) {
    initial <- if(any(h@alpha<0)) h@alpha/2 else .03/ncol(h@beta)
    if(identical(h,models$rejection)) initial <- 2
    compare(h,size,initial,seed)
  }
}
do.call(RNGkind,as.list(saved))
for(h in models) {
  initial <- if(any(h@alpha<0)) h@alpha/2 else .03/ncol(h@beta)
  path <- compare(h,301,initial,913);rng <- .Random.seed
  for(engine in c("r","cpp")) {
    op <- options(emhawkes.engine=engine)
    replay <- logLik(h,path$inter_arrival,path$type,path$mark,lambda_component0=initial,
      N0=2,Nc0=5,infer=TRUE)
    options(op)
    for(field in c("lambda","rambda","lambda_component","rambda_component"))
      stopifnot(isTRUE(all.equal(as.numeric(path[[field]]),as.numeric(replay[[field]]),tolerance=1e-12)))
    stopifnot(identical(rng,.Random.seed))
  }
}
# Automatic initialization and a zero excitation row still consume exp/uniform pairs.
h <- models[["2:2"]]
for(initial in list(NULL,0,matrix(0,2,2))) compare(h,101,initial,511)
if(exists(".Random.seed",.GlobalEnv)) rm(".Random.seed",envir=.GlobalEnv)
hsim(h,1,lambda_component0=0)
stopifnot(!exists(".Random.seed",.GlobalEnv))
# Signed-state rejection after an event cannot restart or discard draws.
h <- new("hspec",mu=.4,alpha=matrix(c(-.2,0),1),beta=matrix(c(1,2),1),type_col_map=list(1:2))
compare(h,5,0,19)
# Match R-level RNG warnings as well as state on reciprocal/bound overflow.
for(level in c(1e-320,1e-300,1e308)) {
  h <- new("hspec",mu=level,alpha=matrix(0,1,2),beta=matrix(c(1,2),1),type_col_map=list(1:2))
  compare(h,3,0,514)
}
h <- new("hspec",mu=1,alpha=matrix(0,1,2),beta=matrix(c(1,2),1),type_col_map=list(1:2))
compare(h,3,1e308,515)
# The existing signed-state tolerance can admit a zero/negative rounded bound.
h <- new("hspec",mu=1e-20,alpha=matrix(0,1,2),beta=matrix(c(1,2),1),type_col_map=list(1:2))
for(level in c(-1e-20,-2e-20)) compare(h,3,matrix(c(level,0),1),516)
# An unequal-rate numeric model actually selects the native wrapper.
.native_calls <- 0L
trace("hawkes_simulate_cpp",where=asNamespace("emhawkes"),
  tracer=quote(.GlobalEnv$.native_calls <- .GlobalEnv$.native_calls+1L),print=FALSE)
hsim(models[["2:2"]],2,lambda_component0=0)
untrace("hawkes_simulate_cpp",where=asNamespace("emhawkes"))
stopifnot(.native_calls==1L)
