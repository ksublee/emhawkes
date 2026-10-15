# Prepare fixtures with tools/benchmark-bid-ask-numeric.R prepare OLD_LIB FIXTURES.
# Rscript --vanilla tools/benchmark-state-gradient.R worker LIB FIXTURES OUT TAG ROUND
# TAG: old, numeric, analytic. Alternate numeric/analytic/old order across rounds.
# Rscript --vanilla tools/benchmark-state-gradient.R compare OUT
args <- commandArgs(TRUE)
if (args[1] == 'compare') {
  out <- args[2]; results <- lapply(list.files(out, pattern='[.]rds$', full.names=TRUE), readRDS)
  stopifnot(length(results) == 9L)
  rows <- lapply(results, function(x) x$timing)
  timing <- do.call(rbind, rows)
  ref <- Filter(function(x) x$tag == 'old' && x$round == 1L, results)[[1]]$fit
  accuracy <- list()
  for (x in results) {
    if (x$tag != 'analytic') for (field in names(ref)) stopifnot(identical(x$fit[[field]], ref[[field]]))
    stopifnot(x$fit$code == 0, abs(x$fit$maximum-ref$maximum) < 2e-6,
      max(abs(x$fit$estimate-ref$estimate)/pmax(1,abs(ref$estimate))) < 2e-5)
    accuracy[[length(accuracy)+1L]] <- data.frame(tag=x$tag,round=x$round,
      maximum_difference=x$fit$maximum-ref$maximum,
      max_parameter_difference=max(abs(x$fit$estimate-ref$estimate)),
      max_relative_parameter_difference=max(abs(x$fit$estimate-ref$estimate)/pmax(1,abs(ref$estimate))),
      code=x$fit$code,iterations=unname(x$fit$iterations))
  }
  summary <- aggregate(seconds~tag,timing,median)
  write.csv(timing,file.path(out,'timings.csv'),row.names=FALSE)
  write.csv(summary,file.path(out,'summary.csv'),row.names=FALSE)
  write.csv(do.call(rbind,accuracy),file.path(out,'accuracy.csv'),row.names=FALSE)
  print(summary,row.names=FALSE);print(do.call(rbind,accuracy),row.names=FALSE);quit()
}
stopifnot(args[1]=='worker')
.libPaths(c(args[2],.libPaths()));library(emhawkes)
stopifnot(identical(normalizePath(dirname(find.package('emhawkes'))),normalizePath(args[2])))
f <- readRDS(file.path(args[3],'fixtures.rds'))$bid_ask_current
a <- f$arguments; out <- args[4]; tag <- args[5]; round <- as.integer(args[6])
dir.create(out,recursive=TRUE,showWarnings=FALSE)
affinity<-Sys.getenv('EMHAWKES_BENCH_AFFINITY','scheduler')
if(affinity!='scheduler') {
 stopifnot(.Platform$OS.type=='windows',grepl('^[1-9][0-9]*$',affinity))
 actual<-system2('powershell.exe',c('-NoProfile','-NonInteractive','-Command',shQuote(sprintf(
  '$p = Get-Process -Id %d; $p.ProcessorAffinity = [IntPtr]%s; $p.ProcessorAffinity.ToInt64()',Sys.getpid(),affinity))),stdout=TRUE)
 stopifnot(length(actual)==1L,trimws(actual)==affinity)
}
options(emhawkes.engine='auto')
if(tag=='analytic') {
 source('tools/state-gradient-cases.R')
 env<-environment(a$object@mu)
 a$gradient<-'analytic';a$jacobian<-bid_ask_jacobians(env$initial_ask_price,
   env$initial_bid_price,env$initial_level,env$initial_mid_price)
}
warm<-a;warm$control<-list(iterlim=2)
for(name in c('inter_arrival','type','mark'))warm[[name]]<-warm[[name]][1:101]
invisible(suppressWarnings(suppressMessages(do.call(hfit,warm))))
gc();elapsed<-unname(system.time(v<-suppressMessages(do.call(hfit,a)))[['elapsed']])
fit<-v[c('estimate','maximum','gradient','hessian','code','iterations')]
timing<-data.frame(tag,round,events=f$events,seconds=elapsed,affinity)
saveRDS(list(tag=tag,round=round,fit=fit,timing=timing),file.path(out,paste0(tag,'-',round,'.rds')))
writeLines(trimws(capture.output(sessionInfo()),which='right'),file.path(out,paste0('session-',tag,'.txt')))
print(timing,row.names=FALSE)
