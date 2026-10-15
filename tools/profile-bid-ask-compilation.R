# Rscript --vanilla tools/profile-bid-ask-compilation.R LIB FIXTURE_DIRECTORY OUTPUT_DIRECTORY
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
.libPaths(c(args[1], .libPaths()))
library(emhawkes)
root <- args[2]
f <- readRDS(file.path(root,'fixtures.rds'))
out <- args[3]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
rows <- list()
for(case in names(f)) {
 ref <- NULL
 for(mode in c('original','compiled_mu','compiled_impact','compiled_both')) {
 a <- f[[case]]$arguments
 if(mode %in% c('compiled_mu','compiled_both')) a$object@mu <- compiler::cmpfun(a$object@mu)
 if(mode %in% c('compiled_impact','compiled_both')) a$object@impact <- compiler::cmpfun(a$object@impact)
 map <- emhawkes:::build_parameter_map(a$object,TRUE)
 obs <- do.call(emhawkes:::normalize_observation,a[intersect(names(a),c('object','inter_arrival','type','mark','N','Nc','N0','Nc0'))])
 objective <- emhawkes:::make_hawkes_objective(a$object,obs,a$lambda_component0,map,derivatives=FALSE)
 p <- map$start
 run <- function() vapply(seq_len(20),function(i) {q <- p; q[1] <- p[1]*(1+1e-7*(i%%2));objective$loglik(q)},0)
 values <- run()
 if(is.null(ref)) ref <- values
 stopifnot(identical(ref,values))
 for(round in 1:3) {
 gc(); time <- system.time(values <- run())[['elapsed']]
 rows[[length(rows)+1L]] <- data.frame(case,mode,round,evaluations=20,seconds=time)
 }
 print(tail(rows,1));flush.console()
 }
}
write.csv(do.call(rbind,rows),file.path(out,'compilation.csv'),row.names=FALSE)
