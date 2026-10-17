# Compare actual initial states and outputs using two separately installed versions.
# Run baseline first, then candidate, using the same results directory:
# Rscript tools/initialization-compare.R baseline /path/to/old/library /path/to/results
# Rscript tools/initialization-compare.R candidate /path/to/new/library /path/to/results
args <- commandArgs(TRUE)
if (length(args) != 3L || !args[1] %in% c('baseline', 'candidate'))
  stop('Supply baseline/candidate, the installed package library, and an existing results directory.')
lab <- normalizePath(args[3], winslash = '/', mustWork = TRUE)
.libPaths(c(normalizePath(args[2], winslash = '/', mustWork = TRUE), .libPaths()))
library(emhawkes)
compiler::enableJIT(0)
quiet <- function(expr) suppressWarnings(suppressMessages(force(expr)))
initial <- function(path) matrix(path$lambda_component[1, ], path$hspec@dimens, byrow=TRUE)
capture <- function(expr) {
  messages <- warnings <- character()
  result <- withCallingHandlers(tryCatch(force(expr),error=identity),
    message=function(m){messages <<- c(messages,conditionMessage(m)); invokeRestart('muffleMessage')},
    warning=function(w){warnings <<- c(warnings,conditionMessage(w)); invokeRestart('muffleWarning')})
  list(value=result,messages=messages,warnings=warnings)
}
if (args[1] == 'baseline') {
  models <- list(
    weak=new('hspec',mu=1,alpha=.05,beta=1),
    moderate=new('hspec',mu=1,alpha=.5,beta=1),
    independent=new('hspec',mu=c(1,1),alpha=diag(.4,2),beta=matrix(1,2,2)),
    cross=new('hspec',mu=c(.7,1),alpha=matrix(c(.12,.03,.04,.1),2),beta=matrix(1.5,2,2)),
    unequal=new('hspec',mu=c(.7,1),alpha=matrix(c(.12,.03,.04,.1),2),beta=matrix(c(.7,1.2,1.7,2.1),2)),
    mapped=new('hspec',mu=c(.7,1),alpha=matrix(seq(.01,.08,.01),2,4),beta=matrix(1.5,2,4),type_col_map=list(c(1,2,3),c(2,4))),
    eta=new('hspec',mu=.8,alpha=.1,beta=1,eta=.03),
    functions=new('hspec',mu=function(param=c(m=.8)) matrix(param[1]),alpha=function(param=c(a=.1)) matrix(param[1]),beta=1),
    observed=new('hspec',mu=function(n,N) matrix(.8+.001*N[max(1,n-1),1]),alpha=.1,beta=1),
    impact=new('hspec',mu=.8,alpha=.1,beta=1,impact=function(lambda_component_n) matrix(.01+.01*lambda_component_n)),
    marked=new('hspec',mu=.8,alpha=.1,beta=1,eta=.03,rmark=function() runif(1,.8,1.2)),
    singular=new('hspec',mu=c(.7,1),alpha=diag(1,2),beta=matrix(1,2,2)))
  states <- lapply(models,function(h) capture(hsim(h,1)))
  fixtures <- list()
  for (name in names(models)) {
    h <- models[[name]]
    if (inherits(states[[name]]$value,'error')) next
    for (seed in c(1,37,804)) for (size in c(1,2,201)) {
      set.seed(seed); path <- quiet(hsim(h,size)); rng <- .Random.seed
      fixtures[[length(fixtures)+1L]] <- list(name=name,seed=seed,size=size,path=path,rng=rng,initial=initial(path))
    }
  }
  fits <- list()
  for (name in c('weak','moderate','cross')) for (events in c(200,2000)) {
    h <- models[[name]]; set.seed(910+events); path <- quiet(hsim(h,events+1))
    p <- emhawkes:::build_parameter_map(h,FALSE)$start
    A <- diag(length(p)); B <- rep(0,length(p))
    if (h@dimens == 1L) { A <- rbind(A,c(0,-1,1)); B <- c(B,0) }
    fit <- capture(hfit(h,path$inter_arrival,type=path$type,optimizer='stats',reduced=FALSE,
      constraints=list(ineqA=A,ineqB=B),control=list(maxit=300)))
    if (inherits(fit$value,'error')) { fits[[paste(name,events)]] <- list(error=conditionMessage(fit$value)); next }
    model <- fitted_hspec(fit$value)
    replay <- quiet(logLik(model,path$inter_arrival,type=path$type,infer=TRUE))
    fits[[paste(name,events)]] <- list(name=name,events=events,truth=h,path=path,
      coefficients=coef(fit$value),maximum=as.numeric(logLik(fit$value)),
      fitted=model,initial=matrix(replay$lambda_component[1,],h@dimens,byrow=TRUE),
      code=fit$value$code,warnings=fit$warnings,arguments=list(constraints=list(ineqA=A,ineqB=B)))
  }
  saveRDS(list(models=models,states=states,fixtures=fixtures,fits=fits),file.path(lab,'baseline.rds'))
  cat('Baseline saved:',length(fixtures),'simulation fixtures;',length(fits),'fit cases\n')
} else {
  baseline <- readRDS(file.path(lab,'baseline.rds'))
  rows <- list(); output_matches <- rng_matches <- reference_matches <- 0L; restore_error <- 0
  for (name in names(baseline$models)) {
    h <- baseline$models[[name]]; out <- capture(hsim(h,1))
    stopifnot(!inherits(out$value,'error'),length(out$messages)==1L)
    prior <- baseline$states[[name]]$value
    old <- if (inherits(prior,'error')) NA_real_ else sum(initial(prior))
    rows[[name]] <- data.frame(model=name,old_component_sum=old,new_component_sum=sum(initial(out$value)),message=out$messages)
  }
  for (fixture in baseline$fixtures) {
    h <- baseline$models[[fixture$name]]
    set.seed(fixture$seed); out <- quiet(hsim(h,fixture$size,lambda_component0=fixture$initial))
    stopifnot(identical(out,fixture$path),identical(.Random.seed,fixture$rng))
    output_matches <- output_matches+1L; rng_matches <- rng_matches+1L
    set.seed(fixture$seed); ref <- quiet(emhawkes:::simulate_hawkes_r(h,fixture$size,lambda_component0=fixture$initial))
    stopifnot(identical(out,ref),identical(.Random.seed,fixture$rng)); reference_matches <- reference_matches+1L
    set.seed(fixture$seed); auto <- quiet(hsim(h,fixture$size)); rng <- .Random.seed
    set.seed(fixture$seed); fixed <- quiet(hsim(h,fixture$size,lambda_component0=initial(auto)))
    stopifnot(identical(auto,fixed),identical(.Random.seed,rng))
    for(engine in c('r','auto')) {
      options(emhawkes.engine=engine)
      replay <- quiet(logLik(h,auto$inter_arrival,type=auto$type,mark=auto$mark,infer=TRUE))
      error <- max(abs(c(auto$lambda-replay$lambda,auto$rambda-replay$rambda)))
      stopifnot(error==0); restore_error <- max(restore_error,error)
    }
  }
  write.csv(do.call(rbind,rows),file.path(lab,'initial-states.csv'),row.names=FALSE)
  fit_rows <- list()
  options(emhawkes.engine='auto')
  for (key in names(baseline$fits)) {
    prior <- baseline$fits[[key]]
    if (!is.null(prior$error)) { cat('Baseline fit failed:',key,prior$error,'\n'); next }
    replay <- quiet(logLik(prior$fitted,prior$path$inter_arrival,type=prior$path$type,lambda_component0=prior$initial))
    stopifnot(identical(as.numeric(replay),prior$maximum))
    fit <- capture(do.call(hfit,c(list(object=prior$truth,inter_arrival=prior$path$inter_arrival,type=prior$path$type,
      optimizer='stats',reduced=FALSE,control=list(maxit=300)),prior$arguments)))
    if (inherits(fit$value,'error')) { cat('Candidate fit failed:',key,conditionMessage(fit$value),'\n'); next }
    new_replay <- quiet(do.call(logLik,c(list(object=fitted_hspec(fit$value),inter_arrival=prior$path$inter_arrival,type=prior$path$type),fit_initial_state(fit$value))))
    stopifnot(identical(as.numeric(new_replay),as.numeric(logLik(fit$value))))
    d <- fit_diagnostics(fit$value)$initialization
    for (parameter in names(prior$coefficients)) fit_rows[[length(fit_rows)+1L]] <- data.frame(model=prior$name,events=prior$events,
      parameter=parameter,old=unname(prior$coefficients[parameter]),new=unname(coef(fit$value)[parameter]),
      old_initial_sum=sum(prior$initial),new_initial_sum=sum(d$components),old_code=prior$code,new_code=fit$value$code)
  }
  write.csv(do.call(rbind,fit_rows),file.path(lab,'fits.csv'),row.names=FALSE)
  cat('Exact old-state output matches:',output_matches,'; RNG matches:',rng_matches,
      '; R reference matches:',reference_matches,'; max intensity restoration error:',restore_error,'\n')
  print(do.call(rbind,rows)[,1:3]); print(do.call(rbind,fit_rows))
}
