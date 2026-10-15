# Rscript --vanilla tools/benchmark-callback-memory.R LIB FIXTURES OUT TAG
# Run old then new, separately from timing workers. Large RDS files stay in OUT.
args<-commandArgs(TRUE)
.libPaths(c(args[1],.libPaths()));library(emhawkes)
stopifnot(identical(normalizePath(dirname(find.package("emhawkes"))),normalizePath(args[1])))
out<-args[3];dir.create(out,recursive=TRUE,showWarnings=FALSE)
affinity<-Sys.getenv("EMHAWKES_BENCH_AFFINITY","scheduler")
if(affinity!="scheduler") {
 stopifnot(.Platform$OS.type=="windows",grepl("^[1-9][0-9]*$",affinity))
 actual<-system2("powershell.exe",c("-NoProfile","-NonInteractive","-Command",shQuote(sprintf(
  "$p = Get-Process -Id %d; $p.ProcessorAffinity = [IntPtr]%s; $p.ProcessorAffinity.ToInt64()",Sys.getpid(),affinity))),stdout=TRUE)
 stopifnot(length(actual)==1L,trimws(actual)==affinity)
}

f<-readRDS(file.path(args[2],"fixtures.rds"))
if(args[4]=="old") {
 set.seed(63001);sim<-hsim(f$bid_ask_current$arguments$object,size=10001L,lambda_component0=f$bid_ask_current$arguments$lambda_component0)
 saveRDS(sim,file.path(out,"data.rds"))
} else sim<-readRDS(file.path(out,"data.rds"))
options(emhawkes.engine="auto");rows<-list()
for(case in names(f))for(n in c(1000L,10000L))for(infer in c(FALSE,TRUE)) {
 a<-list(object=f[[case]]$arguments$object,inter_arrival=sim$inter_arrival[seq_len(n+1L)],type=sim$type[seq_len(n+1L)],mark=sim$mark[seq_len(n+1L)],lambda_component0=f[[case]]$arguments$lambda_component0,infer=infer)
 invoke<-function()suppressMessages(do.call(logLik,a))
 invisible(invoke());gc(reset=TRUE);profile<-tempfile()
 Rprofmem(profile);value<-tryCatch(invoke(),finally=Rprofmem(NULL));heap<-gc()
 path<-file.path(out,paste(case,n,infer,"rds",sep="."))
 if(args[4]=="old")saveRDS(value,path)else stopifnot(identical(value,readRDS(path)))
 con<-file(profile,"r");bytes<-allocations<-0
 repeat {lines<-readLines(con,n=10000L,warn=FALSE);if(!length(lines))break
  sizes<-suppressWarnings(as.numeric(sub(" .*","",substr(lines,1,24))));bytes<-bytes+sum(sizes,na.rm=TRUE);allocations<-allocations+sum(!is.na(sizes))}
 close(con);unlink(profile)
 rows[[length(rows)+1L]]<-data.frame(case,events=n,infer,tag=args[4],allocated_bytes=bytes,allocations,max_R_heap_MiB=sum(heap[,6]),exact=TRUE)
}
write.csv(do.call(rbind,rows),file.path(out,paste0(args[4],".csv")),row.names=FALSE)
