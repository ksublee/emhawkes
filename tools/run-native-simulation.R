# Rscript --vanilla tools/run-native-simulation.R BASE_LIB NEW_LIB RAW OUT
a <- commandArgs(TRUE);stopifnot(length(a)==4L)
libs <- setNames(normalizePath(a[1:2],winslash="/",mustWork=TRUE),c("base","native"))
raw <- a[3];out <- a[4]
dir.create(raw,recursive=TRUE,showWarnings=FALSE);dir.create(out,recursive=TRUE,showWarnings=FALSE)
rscript <- file.path(R.home("bin"),if(.Platform$OS.type=="windows")"Rscript.exe"else"Rscript")
run <- function(script,args,log) {
  status <- system2(rscript,c("--vanilla",shQuote(script),shQuote(args)),stdout=log,stderr=log)
  if(status!=0)stop("Worker failed: ",log)
}
# Cross-version correctness includes callbacks, residuals, extinction and seeds.
for(label in names(libs)) run("tools/check-simulation-copy.R",c(libs[[label]],raw,label),file.path(raw,paste0(label,".checks.log")))
old <- readRDS(file.path(raw,"base.checks.rds"));new <- readRDS(file.path(raw,"native.checks.rds"))
stopifnot(identical(old,new))
write.csv(data.frame(case=names(old),identical=TRUE),file.path(out,"cross-version-checks.csv"),row.names=FALSE)
rm(old,new);gc()
# Sampling is separate from elapsed-time and allocation measurements.
for(label in names(libs)) run("tools/benchmark-native-simulation.R",
  c(libs[[label]],raw,label,"multi_equal",100000,0,"profile"),file.path(raw,paste0(label,".profile.log")))
design <- rbind(expand.grid(case=c("univariate","multivariate","multi_equal","large_equal"),
  events=c(100L,10000L,100000L),stringsAsFactors=FALSE),
  data.frame(case="univariate",events=c(0L,10L)),
  data.frame(case=c("multi_unequal","marked","state_current"),events=10000L))
write.csv(design,file.path(out,"design.csv"),row.names=FALSE)
metrics <- list();counter <- 0L
for(i in seq_len(nrow(design))) {
 z <- design[i,]
 for(round in 1:3) {
   order <- if((i+round)%%2L)c("base","native")else c("native","base")
   values <- list()
   for(label in order) {
     key<-paste(label,z$case,z$events,round,"timing",sep=".")
     run("tools/benchmark-native-simulation.R",c(libs[[label]],raw,label,z$case,z$events,round,"timing"),file.path(raw,paste0(key,".log")))
     values[[label]]<-readRDS(file.path(raw,paste0(key,".rds")))
     counter<-counter+1L;metrics[[counter]]<-values[[label]]$metrics
   }
   stopifnot(identical(values$base$data,values$native$data),identical(values$base$seed,values$native$seed))
   write.csv(do.call(rbind,metrics),file.path(out,"measurements.csv"),row.names=FALSE)
   rm(values);gc()
   cat(z$case,z$events,"round",round,"exact data/RNG match\n");flush.console()
 }
}
# Peak RSS and cumulative allocation are collected in fresh, separate workers.
for(id in c("univariate","multi_equal","large_equal","state_current")) {
 values<-list()
 for(label in names(libs)) {
   key<-paste(label,id,10000,0,"memory",sep=".")
   run("tools/benchmark-native-simulation.R",c(libs[[label]],raw,label,id,10000,0,"memory"),file.path(raw,paste0(key,".log")))
   values[[label]]<-readRDS(file.path(raw,paste0(key,".rds")))
   if(label=="native") worker_session <- values[[label]]$session
   counter<-counter+1L;metrics[[counter]]<-values[[label]]$metrics
 }
 stopifnot(identical(values$base$data,values$native$data),identical(values$base$seed,values$native$seed))
 rm(values);gc()
}
measurements <- do.call(rbind,metrics)
write.csv(measurements,file.path(out,"measurements.csv"),row.names=FALSE)
timings <- subset(measurements,mode=="timing")
summary <- aggregate(cbind(elapsed,cpu)~version+case+events,timings,median)
wide <- merge(subset(summary,version=="base"),subset(summary,version=="native"),by=c("case","events"),suffixes=c("_base","_native"))
wide$speedup <- wide$elapsed_base/wide$elapsed_native
write.csv(wide,file.path(out,"timing-summary.csv"),row.names=FALSE)
write.csv(subset(measurements,mode=="memory"),file.path(out,"memory.csv"),row.names=FALSE)
for(label in names(libs)) {
 file.copy(file.path(raw,paste0(label,".profile-self.csv")),out,overwrite=TRUE)
 file.copy(file.path(raw,paste0(label,".profile-total.csv")),out,overwrite=TRUE)
}
writeLines(c(paste("baseline", "ac880e1f55f6ae18a4c33c8fe32a6928901115df"),
  paste("processor", Sys.getenv("PROCESSOR_IDENTIFIER", Sys.info()[["machine"]])),
  "Three alternating-order rounds; separate processes, sequential workers.",
  "Elapsed uses high-resolution Sys.time(); CPU uses proc.time().",
  "Small requests use repeated complete hsim calls; setup and installation excluded.",
  "Rprofmem is cumulative R allocation; process peak RSS includes loading/warm-up."),file.path(out,"manifest.txt"))
writeLines(sub("[[:blank:]]+$","",worker_session),file.path(out,"session-info.txt"))
cat("All correctness and timing comparisons passed.\n")
