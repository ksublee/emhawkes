# Rscript --vanilla tools/check-state-gradient.R LIB FIXTURES RESULTS
args <- commandArgs(TRUE)
.libPaths(c(args[1],.libPaths()));library(emhawkes)
source('tools/state-gradient-cases.R')
a <- readRDS(file.path(args[2],'fixtures.rds'))$bid_ask_current$arguments
out <- args[3];dir.create(out,recursive=TRUE,showWarnings=FALSE)
env <- environment(a$object@mu)
jac <- bid_ask_jacobians(env$initial_ask_price,env$initial_bid_price,env$initial_level,env$initial_mid_price)
map <- emhawkes:::build_parameter_map(a$object,TRUE)
obs <- do.call(emhawkes:::normalize_observation,a[names(a)!='lambda_component0'])
analytic <- emhawkes:::make_hawkes_objective(a$object,obs,a$lambda_component0,map,jacobian=jac)
numeric <- emhawkes:::make_hawkes_objective(a$object,obs,a$lambda_component0,map,derivatives=FALSE)
options(emhawkes.engine='auto')
fits <- list(); comparisons <- list()
for(start in c('lower','higher')) {
 p <- map$start * if(start=='lower') .9 else 1.1
 current <- a;v <- emhawkes:::restore_parameters(map,p)
 for(name in names(v))slot(current$object,name)<-v[[name]]
 current$control <- list(reltol=1e-12,iterlim=1000)
 fn <- suppressMessages(do.call(hfit,current))
 fa <- suppressMessages(do.call(hfit,c(current,list(gradient='analytic',jacobian=jac))))
 stopifnot(fn$code==0,fa$code==0,abs(fn$maximum-fa$maximum)<1e-6,
  max(abs(coef(fn)-coef(fa))/pmax(1,abs(coef(fn))))<2e-4)
 fits[[start]]<-list(numeric=fn,analytic=fa)
 comparisons[[start]]<-data.frame(start,numeric_loglik=fn$maximum,analytic_loglik=fa$maximum,
  loglik_difference=fa$maximum-fn$maximum,max_parameter_difference=max(abs(coef(fn)-coef(fa))),
  max_relative_parameter_difference=max(abs(coef(fn)-coef(fa))/pmax(1,abs(coef(fn)))),
  numeric_code=fn$code,analytic_code=fa$code)
 print(comparisons[[start]],row.names=FALSE);flush.console()
}
points <- list(start=map$start,nearby=map$start*c(1.01,.99,rep(1.005,length(map$start)-2)),
  numeric_fitted=coef(fits$lower$numeric),analytic_fitted=coef(fits$higher$analytic))
rows <- list()
for(name in names(points)) {
 p <- points[[name]]
 g <- analytic$gradient(p)
 values <- emhawkes:::restore_parameters(map,p);candidate<-a
 for(component in names(values))slot(candidate$object,component)<-values[[component]]
 options(emhawkes.engine='r');reference<-suppressMessages(as.numeric(do.call(logLik,candidate)));options(emhawkes.engine='auto')
 stopifnot(abs(reference-analytic$loglik(p))<1e-9,identical(analytic$loglik(p),numeric$loglik(p)))
 for(step in c(2e-6,1e-6,5e-7)) {
  fd <- vapply(seq_along(p),function(k){up<-down<-p;h<-step*max(1,abs(p[k]));up[k]<-up[k]+h;down[k]<-down[k]-h;
   (numeric$loglik(up)-numeric$loglik(down))/(2*h)},0)
  error <- max(abs(g-fd)/pmax(1,abs(fd)))
  stopifnot(error<2e-5)
  rows[[length(rows)+1L]]<-data.frame(point=name,step,loglik_R=reference,
    loglik_cpp=unname(analytic$loglik(p)),max_absolute_gradient_difference=max(abs(g-fd)),
    max_scaled_gradient_difference=error)
 }
}
write.csv(do.call(rbind,comparisons),file.path(out,'multiple-starts.csv'),row.names=FALSE)
write.csv(do.call(rbind,rows),file.path(out,'gradient-checks.csv'),row.names=FALSE)
saveRDS(fits,file.path(out,'multiple-start-fits.private'))
print(do.call(rbind,rows),row.names=FALSE)
