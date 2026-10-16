# Rscript --vanilla tools/check-simulation-unequal.R NEW_LIB OUT
a <- commandArgs(TRUE);stopifnot(length(a)==2L)
.libPaths(c(a[1],.libPaths()));library(emhawkes)
source("tools/simulation-unequal-cases.R")
r_sim <- getFromNamespace("simulate_hawkes_r","emhawkes")
rows <- draws <- list()
for(id in c("small_unequal","multi_unequal","large_unequal","rejection","signed_unequal")) {
  x <- simulation_profile_case(id)
  set.seed(x$seed);path <- hsim(x$object,1001,lambda_component0=x$initial,N0=2,Nc0=5);rng <- .Random.seed
  .exp_calls <- .uniform_calls <- 0L
  trace("rexp",where=asNamespace("stats"),tracer=quote(.GlobalEnv$.exp_calls <- .GlobalEnv$.exp_calls+1L),print=FALSE)
  trace("runif",where=asNamespace("stats"),tracer=quote(.GlobalEnv$.uniform_calls <- .GlobalEnv$.uniform_calls+1L),print=FALSE)
  set.seed(x$seed);reference <- r_sim(x$object,1001,lambda_component0=x$initial,N0=2,Nc0=5)
  untrace("rexp",where=asNamespace("stats"));untrace("runif",where=asNamespace("stats"))
  stopifnot(identical(path,reference),identical(rng,.Random.seed),.exp_calls==.uniform_calls)
  draws[[id]] <- data.frame(case=id,events=1000L,dimension=x$object@dimens,
    exponential_draws=.exp_calls,uniform_draws=.uniform_calls,
    rejected_candidates=.exp_calls-1000L*x$object@dimens)
  replay <- list()
  for(engine in c("r","cpp")) {
    op <- options(emhawkes.engine=engine)
    replay[[engine]] <- logLik(x$object,path$inter_arrival,path$type,path$mark,
      lambda_component0=x$initial,N0=2,Nc0=5,infer=TRUE)
    options(op)
    for(field in c("lambda","lambda_component","rambda","rambda_component")) {
      error <- max(abs(path[[field]]-replay[[engine]][[field]])/pmax(1,abs(path[[field]])))
      stopifnot(error<=1e-12)
      rows[[length(rows)+1L]] <- data.frame(case=id,engine=engine,field=field,scaled_error=error)
    }
    stopifnot(identical(rng,.Random.seed))
  }
  error <- abs(replay$r$loglikelihood-replay$cpp$loglikelihood)/1000
  stopifnot(error<=1e-12)
  rows[[length(rows)+1L]] <- data.frame(case=id,engine="r_cpp",field="loglikelihood_per_event",scaled_error=error)
}
write.csv(do.call(rbind,rows),file.path(a[2],"replay.csv"),row.names=FALSE)
write.csv(do.call(rbind,draws),file.path(a[2],"draws.csv"),row.names=FALSE)
cat("Five thinning fixtures match exactly, including rejected draws and independent inference.\n")
