# Rscript --vanilla tools/benchmark-callback-overhead.R worker LIB FIXTURES OUT TAG ROUND [CASE]
# Run isolated workers sequentially, alternating old/new order across three rounds.
# Rscript --vanilla tools/benchmark-callback-overhead.R compare OUT
args <- commandArgs(TRUE)
if(args[1]=="compare") {
 out <- args[2]; files <- list.files(out,pattern="[.]rds$",full.names=TRUE)
 data <- lapply(files,readRDS); rows <- checks <- list()
 stopifnot(length(data)==12L)
 for(case in c("bid_ask","bid_ask_current"))for(tag in c("old","new")) {
  rounds<-vapply(Filter(function(x)x$case==case && x$tag==tag,data),`[[`,0L,"round")
  stopifnot(identical(sort(rounds),1:3))
 }
 for(case in c("bid_ask","bid_ask_current")) {
  ref <- Filter(function(x)x$case==case && x$tag=="old" && x$round==1L,data)[[1]]
  for(x in Filter(function(x)x$case==case,data)) {
   for(field in names(ref$values)) {
    a<-ref$values[[field]];b<-x$values[[field]]
    stopifnot(identical(a,b))
    for(item in if(field=="fit")names(a) else field) {
     left<-if(field=="fit")a[[item]] else a
     right<-if(field=="fit")b[[item]] else b
     checks[[length(checks)+1L]] <- data.frame(case,tag=x$tag,round=x$round,
      field=if(field=="fit")paste0("fit.",item) else field,exact=TRUE,
      max_abs_difference=max(abs(unlist(left,use.names=FALSE)-unlist(right,use.names=FALSE))))
    }
   }
   rows[[length(rows)+1L]]<-x$timings
  }
 }
 timings<-do.call(rbind,rows);summary<-aggregate(seconds~case+mode+tag,timings,median)
 write.csv(timings,file.path(out,"timings.csv"),row.names=FALSE)
 write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE)
 write.csv(do.call(rbind,checks),file.path(out,"accuracy.csv"),row.names=FALSE)
 print(summary,row.names=FALSE);quit()
}
stopifnot(args[1]=="worker")
.libPaths(c(args[2],.libPaths()));library(emhawkes)
stopifnot(identical(normalizePath(dirname(find.package("emhawkes"))),normalizePath(args[2])))
f<-readRDS(file.path(args[3],"fixtures.rds"));out<-args[4];tag<-args[5];round<-as.integer(args[6])
dir.create(out,recursive=TRUE,showWarnings=FALSE)
affinity<-Sys.getenv("EMHAWKES_BENCH_AFFINITY","scheduler")
if(affinity!="scheduler") {
 stopifnot(.Platform$OS.type=="windows",grepl("^[1-9][0-9]*$",affinity))
 actual<-system2("powershell.exe",c("-NoProfile","-NonInteractive","-Command",shQuote(sprintf(
  "$p = Get-Process -Id %d; $p.ProcessorAffinity = [IntPtr]%s; $p.ProcessorAffinity.ToInt64()",Sys.getpid(),affinity))),stdout=TRUE)
 stopifnot(length(actual)==1L,trimws(actual)==affinity)
}
if(length(args)>=7L)f<-f[args[7]]
options(emhawkes.engine="auto")
for(case in names(f)) {
 a<-f[[case]]$arguments
 obs<-do.call(emhawkes:::normalize_observation,a[intersect(names(a),c("object","inter_arrival","type","mark","N","Nc","N0","Nc0"))])
 map<-emhawkes:::build_parameter_map(a$object,TRUE)
 objective<-emhawkes:::make_hawkes_objective(a$object,obs,a$lambda_component0,map,derivatives=FALSE)
 p<-map$start
 public<-a[intersect(names(a),c("object","inter_arrival","type","mark","N","Nc","N0","Nc0","lambda_component0"))]
 values<-timings<-list()
 # Warm callback bodies and optimizer setup before recording elapsed times.
 invisible(suppressMessages(do.call(logLik,public)))
 warm<-a;warm$control<-list(iterlim=2)
 for(name in intersect(names(warm),c("inter_arrival","type","mark")))warm[[name]]<-warm[[name]][seq_len(min(101L,length(warm[[name]])))]
 for(name in intersect(names(warm),c("N","Nc")))warm[[name]]<-warm[[name]][seq_len(min(101L,nrow(warm[[name]]))),,drop=FALSE]
 invisible(suppressWarnings(suppressMessages(do.call(hfit,warm))))
 for(mode in c("objective","loglik","inference","fit")) {
  batch<-if(mode=="fit")1L else 20L
  run<-switch(mode,
   objective=function()vapply(seq_len(batch),function(i){q<-p;q[1]<-p[1]*(1+1e-7*(i%%2));objective$loglik(q)},0),
   loglik=function(){for(i in seq_len(batch))v<-suppressMessages(do.call(logLik,public));v},
   inference=function(){for(i in seq_len(batch))v<-suppressMessages(do.call(logLik,c(public,list(infer=TRUE))));v},
   fit=function(){v<-suppressWarnings(suppressMessages(do.call(hfit,a)));v[c("estimate","maximum","gradient","hessian","code","iterations")]})
  if(mode!="fit")invisible(run())
  gc();elapsed<-unname(system.time(value<-run())[["elapsed"]])
  values[[mode]]<-value
  timings[[mode]]<-data.frame(case,events=f[[case]]$events,mode,tag,round,batch,seconds=elapsed/batch,affinity)
  cat(tag,round,case,mode,elapsed/batch,"seconds\n");flush.console()
 }
 saveRDS(list(case=case,tag=tag,round=round,values=values,timings=do.call(rbind,timings)),
  file.path(out,paste(tag,round,case,"rds",sep=".")))
 writeLines(trimws(capture.output(sessionInfo()),which="right"),file.path(out,paste0("session-",tag,".txt")))
}
