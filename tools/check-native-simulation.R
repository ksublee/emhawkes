# Rscript --vanilla tools/check-native-simulation.R NEW_LIB OUT
a <- commandArgs(TRUE);.libPaths(c(a[1],.libPaths()));library(emhawkes)
source("tools/simulation-profile-cases.R")
r_sim <- getFromNamespace("simulate_hawkes_r","emhawkes")
dir.create(a[2],recursive=TRUE,showWarnings=FALSE)
rows <- list()
for(id in c("univariate","multivariate","multi_equal","multi_unequal","large_equal",
  "large_unequal","marked","parameter_mu","state_current","state_dots","state_observed","signed")) {
  x <- simulation_profile_case(id);set.seed(x$seed)
  path <- hsim(x$object,size=1001,lambda_component0=x$initial);rng <- .Random.seed
  set.seed(x$seed);reference <- r_sim(x$object,size=1001,lambda_component0=x$initial)
  stopifnot(identical(rng,.Random.seed),
    identical(path[setdiff(names(path),"hspec")],reference[setdiff(names(reference),"hspec")]))
  replay <- list()
  for(engine in c("r","cpp")) {
    op <- options(emhawkes.engine=engine)
    replay[[engine]] <- suppressMessages(logLik(x$object,path$inter_arrival,path$type,path$mark,
      lambda_component0=x$initial,infer=TRUE));options(op)
    stopifnot(identical(rng,.Random.seed))
    for(field in c("lambda","lambda_component","rambda","rambda_component")) {
      error <- max(abs(path[[field]]-replay[[engine]][[field]])/pmax(1,abs(path[[field]])))
      stopifnot(error<=1e-12)
      rows[[length(rows)+1L]] <- data.frame(case=id,engine,field,scaled_error=error)
    }
  }
  error <- abs(replay$r$loglikelihood-replay$cpp$loglikelihood)/1000
  stopifnot(error<=1e-12)
  rows[[length(rows)+1L]] <- data.frame(case=id,engine="r_cpp",field="loglikelihood_per_event",scaled_error=error)
}
write.csv(do.call(rbind,rows),file.path(a[2],"replay.csv"),row.names=FALSE)
cat("All 12 fixtures, including fallback callbacks, pass independent inference replay.\n")
