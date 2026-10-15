# Rscript --vanilla tools/check-gradient-scaling-fits.R LIB OUT
args<-commandArgs(TRUE);.libPaths(c(args[1],.libPaths()));library(emhawkes)
stopifnot(identical(normalizePath(dirname(find.package('emhawkes'))),normalizePath(args[1])))
out<-args[2];grid<-read.csv(file.path(out,'design.csv'),stringsAsFactors=FALSE)
rows<-list()
for(id in grid$id[grid$fit]) {
  x<-readRDS(file.path(out,'fixtures',paste0(id,'.rds')));a<-x$arguments
  saved<-readRDS(file.path(out,'samples',paste(id,'fit','numeric',1,'rds',sep='.')))
  events<-saved$meta$events_used
  for(name in c('inter_arrival','type','mark'))a[[name]]<-a[[name]][seq_len(events+1L)]
  x$meta$events_used<-events
  obs<-do.call(emhawkes:::normalize_observation,a[intersect(names(a),c('object','inter_arrival','type','mark'))])
  map<-emhawkes:::build_parameter_map(a$object,TRUE)
  obj<-emhawkes:::make_hawkes_objective(a$object,obs,a$lambda_component0,map,jacobian=x$jacobian)
  for(method in c('numeric','analytic')) {
    result<-readRDS(file.path(out,'samples',paste(id,'fit',method,1,'rds',sep='.')))$value
    candidate<-a$object;values<-emhawkes:::restore_parameters(map,result$estimate)
    for(name in names(values))slot(candidate,name)<-values[[name]]
    options(emhawkes.engine='r')
    reference<-suppressMessages(as.numeric(logLik(candidate,a$inter_arrival,a$type,a$mark,lambda_component0=a$lambda_component0)))
    options(emhawkes.engine='auto')
    score<-obj$gradient(result$estimate)
    direction<-cos(seq_along(result$estimate));direction<-direction/sqrt(sum(direction^2))
    directional_fd<-(obj$loglik(result$estimate+2e-5*direction)-
      obj$loglik(result$estimate-2e-5*direction))/4e-5
    directional_error<-abs(sum(score*direction)-directional_fd)/max(1,abs(directional_fd))
    stopifnot(directional_error<5e-5)
    H<-result$hessian;info<--(H+t(H))/2
    scale<-sqrt(abs(diag(info)))
    eigenvalues<-if(all(is.finite(info))&&all(scale>0))eigen(info/outer(scale,scale),symmetric=TRUE,only.values=TRUE)$values else NA_real_
    smallest<-min(eigenvalues);condition<-if(is.finite(smallest)&&smallest>0)max(eigenvalues)/smallest else Inf
    stopifnot(abs(reference-result$maximum)<1e-7)
    rows[[length(rows)+1L]]<-cbind(x$meta,data.frame(method,loglik_R=reference,loglik_fitted=result$maximum,
      loglik_abs_difference=abs(reference-result$maximum),analytic_gradient_infinity=max(abs(score)),
      analytic_gradient_per_event=max(abs(score))/events,
      fitted_directional_scaled_difference=directional_error,
      stored_gradient_difference=max(abs(result$gradient-score)),
      scaled_information_smallest_eigenvalue=smallest,scaled_information_condition=condition))
  }
  cat('Verified fitted likelihood',id,'\n');flush.console()
  write.csv(do.call(rbind,rows),file.path(out,'fitted-reference-checks.csv'),row.names=FALSE)
}
