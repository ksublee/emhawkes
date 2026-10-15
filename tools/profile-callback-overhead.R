# Rscript --vanilla tools/profile-callback-overhead.R LIB FIXTURES OUT [CASE]
args <- commandArgs(TRUE)
.libPaths(c(args[1], .libPaths()))
library(emhawkes)
dir.create(args[3], recursive=TRUE, showWarnings=FALSE)
f <- readRDS(file.path(args[2], "fixtures.rds"))
affinity<-Sys.getenv("EMHAWKES_BENCH_AFFINITY","scheduler")
if(affinity!="scheduler") {
 stopifnot(.Platform$OS.type=="windows",grepl("^[1-9][0-9]*$",affinity))
 actual<-system2("powershell.exe",c("-NoProfile","-NonInteractive","-Command",shQuote(sprintf(
  "$p = Get-Process -Id %d; $p.ProcessorAffinity = [IntPtr]%s; $p.ProcessorAffinity.ToInt64()",Sys.getpid(),affinity))),stdout=TRUE)
 stopifnot(length(actual)==1L,trimws(actual)==affinity)
}
if(length(args)>=4L)f<-f[args[4]]
rows <- list()
for(case in names(f)) {
 a <- f[[case]]$arguments
 obs <- do.call(emhawkes:::normalize_observation, a[intersect(names(a),c("object","inter_arrival","type","mark","N","Nc","N0","Nc0"))])
 map <- emhawkes:::build_parameter_map(a$object,TRUE)
 objective <- emhawkes:::make_hawkes_objective(a$object,obs,a$lambda_component0,map,derivatives=FALSE)
 p <- map$start
 public <- a[intersect(names(a),c("object","inter_arrival","type","mark","N","Nc","N0","Nc0","lambda_component0"))]
 options(emhawkes.engine="auto")
 for(mode in c("objective", "loglik", "inference", "fit")) {
  run <- switch(mode,
   objective=function() for(i in 1:100) {q <- p; q[1] <- p[1]*(1+1e-7*(i%%2)); objective$loglik(q)},
   loglik=function() for(i in 1:100) suppressMessages(do.call(logLik,public)),
   inference=function() for(i in 1:100) suppressMessages(do.call(logLik,c(public,list(infer=TRUE)))),
   fit=function() suppressWarnings(suppressMessages(do.call(hfit,a))))
  if(mode!="fit") invisible(run())
  gc(); path <- file.path(args[3],paste0(case,"-",mode,".Rprof"))
  Rprof(path,interval=.001)
  elapsed <- tryCatch(unname(system.time(value <- run())[["elapsed"]]),finally=Rprof(NULL))
  s <- summaryRprof(path)
  write.csv(s$by.self, file.path(args[3],paste0(case,"-",mode,"-self.csv")))
  write.csv(s$by.total, file.path(args[3],paste0(case,"-",mode,"-total.csv")))
  rows[[length(rows)+1L]] <- data.frame(case,mode,evaluations=if(mode=="fit")1L else 100L,seconds=elapsed)
  write.csv(do.call(rbind,rows),file.path(args[3],"profile-timings.csv"),row.names=FALSE)
  if(mode=="fit") saveRDS(value[c("estimate","maximum","gradient","hessian","code","iterations")],file.path(args[3],paste0(case,"-fit.rds")))
  cat(case,mode,elapsed,"seconds\n");print(head(s$by.self,8));flush.console()
 }
}
