# Rscript --vanilla tools/check-gradient-scaling-hessians.R LIB OUT
# Diagnose default-Hessian failures without changing hfit() or its defaults.
args<-commandArgs(TRUE);.libPaths(c(args[1],.libPaths()));library(emhawkes)
stopifnot(identical(normalizePath(dirname(find.package('emhawkes'))),normalizePath(args[1])))
out<-args[2];options(emhawkes.engine='auto')
checks<-read.csv(file.path(out,'fitted-reference-checks.csv'),stringsAsFactors=FALSE)
ids<-checks$id[checks$method=='numeric' & checks$scaled_information_smallest_eigenvalue<=0]
rows<-list()
for(id in ids) {
  x<-readRDS(file.path(out,'fixtures',paste0(id,'.rds')))
  numeric_fit<-readRDS(file.path(out,'samples',paste(id,'fit','numeric',1,'rds',sep='.')))
  a<-x$arguments;events<-numeric_fit$meta$events_used
  for(name in c('inter_arrival','type','mark'))a[[name]]<-a[[name]][seq_len(events+1L)]
  obs<-do.call(emhawkes:::normalize_observation,a[intersect(names(a),c('object','inter_arrival','type','mark'))])
  map<-emhawkes:::build_parameter_map(a$object,TRUE)
  objective<-emhawkes:::make_hawkes_objective(a$object,obs,a$lambda_component0,map,jacobian=x$jacobian)
  theta<-numeric_fit$value$estimate;previous<-NULL
  for(step in c(1e-4,2e-5)) {
    H<-vapply(seq_along(theta),function(i) {
      plus<-minus<-theta;h<-step*max(1,abs(theta[i]))
      plus[i]<-plus[i]+h;minus[i]<-minus[i]-h
      (objective$gradient(plus)-objective$gradient(minus))/(2*h)
    },numeric(length(theta)))
    info<--(H+t(H))/2;scale<-sqrt(abs(diag(info)))
    ev<-eigen(info/outer(scale,scale),symmetric=TRUE,only.values=TRUE)$values
    difference<-if(is.null(previous))NA_real_ else max(abs(info-previous)/pmax(1,abs(previous)))
    stopifnot(min(ev)>0,is.na(difference)||difference<1e-4)
    rows[[length(rows)+1L]]<-data.frame(id,events_used=events,parameters=length(theta),step,
      scaled_information_smallest_eigenvalue=min(ev),scaled_information_condition=max(ev)/min(ev),
      max_scaled_change_from_previous_step=difference,
      max_scaled_difference_from_default=max(abs(info+numeric_fit$value$hessian)/pmax(1,abs(info))))
    previous<-info
  }
  cat('Checked score-based Hessian',id,'\n');flush.console()
  write.csv(do.call(rbind,rows),file.path(out,'hessian-cross-checks.csv'),row.names=FALSE)
}
