# Rscript --vanilla tools/profile-bid-ask-bound.R LIB FIXTURE_DIRECTORY OUTPUT_DIRECTORY
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
.libPaths(c(args[1],.libPaths()))
library(emhawkes)
f <- readRDS(file.path(args[2],'fixtures.rds'))
ns <- asNamespace('emhawkes'); original <- get('bind_component',ns)
compiled <- function(...) compiler::cmpfun(original(...))
rows <- list()
for(case in names(f)) {
 a <- f[[case]]$arguments; map <- emhawkes:::build_parameter_map(a$object,TRUE)
 obs <- do.call(emhawkes:::normalize_observation,a[intersect(names(a),c('object','inter_arrival','type','mark','N','Nc','N0','Nc0'))])
 reference <- NULL
 for(round in 1:3) for(mode in if(round%%2) c('original','compile_bound') else c('compile_bound','original')) {
 unlockBinding('bind_component',ns);assign('bind_component',if(mode=='original') original else compiled,ns);lockBinding('bind_component',ns)
 objective <- emhawkes:::make_hawkes_objective(a$object,obs,a$lambda_component0,map,derivatives=FALSE)
 p <- map$start
 run <- function() vapply(1:100,function(i) {q<-p;q[1]<-p[1]*(1+1e-7*(i%%2));objective$loglik(q)},0)
 gc();tm<-system.time(v<-run())[['elapsed']]
 if(is.null(reference)) reference<-v
 stopifnot(identical(reference,v))
 rows[[length(rows)+1]]<-data.frame(case,mode,round,evaluations=100,seconds=tm)
 print(tail(rows,1));flush.console()
 }
}
write.csv(do.call(rbind,rows),file.path(args[3],'compile-bound.csv'),row.names=FALSE)
