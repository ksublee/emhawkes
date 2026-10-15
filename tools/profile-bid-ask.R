# Rscript --vanilla tools/profile-bid-ask.R LIB FIXTURE_DIRECTORY OUTPUT_DIRECTORY
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
.libPaths(c(args[1], .libPaths()))
library(emhawkes)
root <- args[2]
fixtures <- readRDS(file.path(root,'fixtures.rds'))
out <- args[3]
dir.create(out,showWarnings=FALSE,recursive=TRUE)
rows <- list()
for (case in names(fixtures)) {
 a <- fixtures[[case]]$arguments
 map <- emhawkes:::build_parameter_map(a$object,TRUE)
 obs <- do.call(emhawkes:::normalize_observation,a[intersect(names(a),c('object','inter_arrival','type','mark','N','Nc','N0','Nc0'))])
 objective <- emhawkes:::make_hawkes_objective(a$object,obs,a$lambda_component0,map,derivatives=FALSE)
 p <- map$start
 run <- function(reps) for(i in seq_len(reps)) { q <- p; q[1] <- p[1]*(1+1e-7*(i%%2)); objective$loglik(q) }
 run(4)
 path <- file.path(out,paste0(case,'.Rprof'))
 Rprof(path,interval=0.001)
 elapsed <- system.time(run(500))[['elapsed']]
 Rprof(NULL)
 s <- summaryRprof(path)
 write.csv(s$by.self,file.path(out,paste0(case,'-self.csv')))
 write.csv(s$by.total,file.path(out,paste0(case,'-total.csv')))
 rows[[case]] <- data.frame(case=case,evaluations=500,seconds=elapsed,milliseconds=elapsed*2)
 print(rows[[case]]); print(head(s$by.self,15)); print(head(s$by.total,15))
}
write.csv(do.call(rbind,rows),file.path(out,'timings.csv'),row.names=FALSE)
