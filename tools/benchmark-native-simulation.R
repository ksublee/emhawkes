# Worker: LIB RAW LABEL CASE EVENTS ROUND MODE [CASES_FILE]
# MODE is timing, memory, rss, or profile. Never run these workers concurrently.
a <- commandArgs(TRUE); stopifnot(length(a) %in% c(7L,8L))
.libPaths(c(a[1],.libPaths())); library(emhawkes)
stopifnot(identical(normalizePath(dirname(find.package("emhawkes"))),normalizePath(a[1])))
source(if(length(a)==8L) a[8] else "tools/simulation-profile-cases.R")
x <- simulation_profile_case(a[4]); events <- as.integer(a[5]); round <- as.integer(a[6]); mode <- a[7]
dir.create(a[2],recursive=TRUE,showWarnings=FALSE)
RNGkind("Mersenne-Twister","Inversion","Rejection")
run <- function(n=events) hsim(x$object,size=n+1L,lambda_component0=x$initial)
for(i in 1:3) {set.seed(x$seed);invisible(run(100L))}
gc();set.seed(x$seed)
iterations <- if(mode=="timing") if(events<=10)1000L else if(events<=100)100L else if(events<=10000)5L else 1L else 1L
allocated <- NA_real_; elapsed <- cpu <- NA_real_
if(mode=="timing") {
  cpu_start <- proc.time()
  # proc.time() elapsed is too coarse for the fastest Windows native calls.
  # Sys.time() uses the high-resolution wall clock; CPU retains proc.time().
  start <- Sys.time()
  for(i in seq_len(iterations)) path <- run()
  elapsed <- as.numeric(difftime(Sys.time(),start,units="secs"))/iterations
  timer <- proc.time()-cpu_start
  cpu <- unname(timer["user.self"]+timer["sys.self"])/iterations
} else if(mode=="memory") {
  file <- tempfile();Rprofmem(file)
  path <- tryCatch(run(),finally=Rprofmem(NULL))
  lines <- readLines(file,warn=FALSE);unlink(file)
  bytes <- suppressWarnings(as.numeric(sub(" .*","",lines)))
  allocated <- sum(bytes[is.finite(bytes)])
} else if(mode=="rss") {
  path <- run()
} else if(mode=="profile") {
  file <- file.path(a[2],paste(a[3],a[4],"Rprof",sep="."))
  Rprof(file,interval=.001,line.profiling=TRUE,gc.profiling=TRUE)
  path <- tryCatch(run(),finally=Rprof(NULL))
  profile <- summaryRprof(file)
  write.csv(profile$by.self,file.path(a[2],paste0(a[3],".profile-self.csv")))
  write.csv(profile$by.total,file.path(a[2],paste0(a[3],".profile-total.csv")))
} else stop("Unknown worker mode")
rng <- .Random.seed
# Query after measurement; the process peak includes package loading and warm-up.
peak <- if(.Platform$OS.type=="windows" && requireNamespace("ps",quietly=TRUE))
  unname(ps::ps_memory_info()["peak_wset"]) else NA_real_
result <- data.frame(version=a[3],case=a[4],events,round,mode,iterations,
  elapsed,cpu,allocated_bytes=allocated,peak_process_rss_bytes=peak,
  output_bytes=as.numeric(object.size(path)))
key <- paste(a[3],a[4],events,round,mode,sep=".")
saveRDS(list(metrics=result,data=path[setdiff(names(path),"hspec")],seed=rng,
  session=capture.output(sessionInfo())),file.path(a[2],paste0(key,".rds")),compress=FALSE)
cat(key,elapsed,"seconds per complete hsim call\n")
