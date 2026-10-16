# Rscript --vanilla tools/run-simulation-unequal.R BASE_LIB NEW_LIB RAW OUT [resume]
# Run alone: compilation, tests and other heavy work must finish first.
a <- commandArgs(TRUE);stopifnot(length(a) %in% c(4L,5L))
resume <- length(a)==5L && identical(a[5],"resume")
libs <- setNames(normalizePath(a[1:2],winslash="/",mustWork=TRUE),c("base","native"))
raw <- a[3];out <- a[4]
dir.create(raw,recursive=TRUE,showWarnings=FALSE);dir.create(out,recursive=TRUE,showWarnings=FALSE)
rscript <- file.path(R.home("bin"),if(.Platform$OS.type=="windows")"Rscript.exe"else"Rscript")
run <- function(script,args,log) {
  status <- system2(rscript,c("--vanilla",shQuote(script),shQuote(args)),stdout=log,stderr=log)
  if(status!=0) stop("Worker failed: ",log)
}
for(label in names(libs)) run("tools/check-simulation-copy.R",c(libs[[label]],raw,label),file.path(raw,paste0(label,".checks.log")))
old <- readRDS(file.path(raw,"base.checks.rds"));new <- readRDS(file.path(raw,"native.checks.rds"))
stopifnot(identical(old,new))
write.csv(data.frame(case=names(old),identical=TRUE),file.path(out,"cross-version-checks.csv"),row.names=FALSE)
rm(old,new);gc()
worker <- function(label,id,events,round,mode) {
  key <- paste(label,id,events,round,mode,sep=".")
  run("tools/benchmark-native-simulation.R",c(libs[[label]],raw,label,id,events,round,mode,
    "tools/simulation-unequal-cases.R"),file.path(raw,paste0(key,".log")))
  readRDS(file.path(raw,paste0(key,".rds")))
}
if(!resume) for(label in names(libs)) worker(label,"multi_unequal",100000L,0L,"profile")
design <- rbind(data.frame(case="small_unequal",events=c(0L,10L,100L,10000L,100000L)),
  data.frame(case="multi_unequal",events=c(100L,10000L,100000L)),
  data.frame(case="large_unequal",events=c(10000L,100000L)),
  data.frame(case=c("rejection","signed_unequal","multi_equal","state_current"),events=10000L),
  data.frame(case="univariate",events=c(0L,10L,100L)))
write.csv(design,file.path(out,"design.csv"),row.names=FALSE)
checkpoint <- file.path(out,"measurements.csv")
previous <- if(resume && file.exists(checkpoint)) read.csv(checkpoint,stringsAsFactors=FALSE) else NULL
metrics <- if(is.null(previous)) list() else lapply(seq_len(nrow(previous)),function(i) previous[i,,drop=FALSE])
completed <- function(id,events,round,mode) {
  if(is.null(previous)) return(FALSE)
  p <- previous[previous$case==id & previous$events==events & previous$round==round & previous$mode==mode,,drop=FALSE]
  if(nrow(p)!=2L || !setequal(p$version,names(libs))) return(FALSE)
  values <- lapply(names(libs),function(label) readRDS(file.path(raw,
    paste0(paste(label,id,events,round,mode,sep="."),".rds"))))
  stopifnot(identical(values[[1]]$data,values[[2]]$data),identical(values[[1]]$seed,values[[2]]$seed))
  TRUE
}
pair <- function(id,events,round,mode,order) {
  values <- list()
  for(label in order) {
    values[[label]] <- worker(label,id,events,round,mode)
    metrics[[length(metrics)+1L]] <<- values[[label]]$metrics
  }
  stopifnot(identical(values$base$data,values$native$data),identical(values$base$seed,values$native$seed))
  write.csv(do.call(rbind,metrics),file.path(out,"measurements.csv"),row.names=FALSE)
  cat(id,events,round,mode,"exact output/RNG match\n");flush.console()
  invisible(NULL)
}
for(i in seq_len(nrow(design))) for(round in 1:3) {
  if(completed(design$case[i],design$events[i],round,"timing")) next
  order <- if((i+round)%%2L)c("base","native")else c("native","base")
  pair(design$case[i],design$events[i],round,"timing",order)
  gc()
}
for(mode in c("memory","rss")) for(id in c("small_unequal","multi_unequal","large_unequal","rejection","signed_unequal","state_current")) {
  if(completed(id,10000L,0L,mode)) next
  pair(id,10000L,0L,mode,names(libs));gc()
}
measurements <- do.call(rbind,metrics)
summary <- aggregate(cbind(elapsed,cpu)~version+case+events,subset(measurements,mode=="timing"),median)
wide <- merge(subset(summary,version=="base"),subset(summary,version=="native"),
  by=c("case","events"),suffixes=c("_base","_native"))
wide$speedup <- wide$elapsed_base/wide$elapsed_native
write.csv(wide,file.path(out,"timing-summary.csv"),row.names=FALSE)
write.csv(subset(measurements,mode %in% c("memory","rss")),file.path(out,"memory.csv"),row.names=FALSE)
for(label in names(libs)) for(kind in c("self","total"))
  file.copy(file.path(raw,paste0(label,".profile-",kind,".csv")),out,overwrite=TRUE)
session <- readRDS(file.path(raw,"native.small_unequal.10000.0.memory.rds"))$session
writeLines(sub("[[:blank:]]+$","",session),file.path(out,"session-info.txt"))
writeLines(c("baseline 53a5f83483146f8efb5363ac62a871c45a6cc7f6",
  paste("processor",Sys.getenv("PROCESSOR_IDENTIFIER",Sys.info()[["machine"]])),
  "Three alternating-order rounds, sequential fresh processes, three 100-event warm-ups.",
  "High-resolution wall time per complete hsim call, setup/installation excluded.",
  "Rprofmem: cumulative R allocation; use separate untraced rss workers for peak working set.",
  "rss peak working set includes process lifetime, package loading and warm-up.",
  "Profile percentages use samples; Windows sample seconds are not elapsed seconds."),file.path(out,"manifest.txt"))
cat("All cross-version comparisons passed.\n")
