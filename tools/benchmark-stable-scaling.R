# Rscript --vanilla tools/benchmark-stable-scaling.R run|compare|reference LIB FIXTURES OUT
# A worker fits one case/method in a fresh process: worker LIB FIXTURES OUT ID METHOD
args <- commandArgs(TRUE)
action <- args[1]; lib <- args[2]; fixtures <- args[3]; out <- args[4]
.libPaths(c(lib,.libPaths())); library(emhawkes)
stopifnot(identical(normalizePath(dirname(find.package('emhawkes'))),normalizePath(lib)))
options(emhawkes.engine='auto')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
checks <- read.csv(file.path(fixtures,'fitted-reference-checks.csv'),stringsAsFactors=FALSE)
design <- checks[checks$method=='numeric' & checks$scaled_information_smallest_eigenvalue<=0,]
design$parameters <- 3L*design$g
design <- design[order(design$parameters,design$feedback,design$observed),]
stopifnot(nrow(design)==5L,all(design$events_used==2500L))
sample_path <- function(id,method)file.path(out,'samples',paste0(id,'.',method,'.rds'))
load_case <- function(id) {
  x<-readRDS(file.path(fixtures,'fixtures',paste0(id,'.rds')))
  events<-design$events_used[match(id,design$id)]
  for(name in c('inter_arrival','type','mark'))x$arguments[[name]]<-x$arguments[[name]][seq_len(events+1L)]
  x$meta$events_used<-events
  x
}
fields <- c('estimate','maximum','gradient','hessian','code','iterations')
information <- function(H) {
  info<--(H+t(H))/2;scale<-sqrt(abs(diag(info)))
  ev<-if(all(is.finite(info))&&all(scale>0))
    eigen(info/outer(scale,scale),symmetric=TRUE,only.values=TRUE)$values else NA_real_
  list(minimum=min(ev),condition=if(all(is.finite(ev))&&min(ev)>0)max(ev)/min(ev) else Inf)
}
if(action=='run') {
  dir.create(file.path(out,'samples'),showWarnings=FALSE)
  write.csv(design,file.path(out,'design.csv'),row.names=FALSE)
  writeLines(c(paste('repository_commit',system2('git',c('rev-parse','HEAD'),stdout=TRUE)),
    'runtime_source 4d6ea845cbce711b3b778bde2f0decfe07a75138',
    paste('R',getRversion()),paste('maxLik',packageVersion('maxLik')),paste('Rcpp',packageVersion('Rcpp'))),file.path(out,'manifest.txt'))
  writeLines(trimws(capture.output(sessionInfo()),which='right'),file.path(out,'session-info.txt'))
  rs<-file.path(R.home('bin'),if(.Platform$OS.type=='windows')'Rscript.exe' else 'Rscript')
  for(i in seq_len(nrow(design)))for(method in if(i%%2)c('maxLik','stable') else c('stable','maxLik')) {
    id<-design$id[i]
    if(file.exists(sample_path(id,method)))next
    cat('Start',id,method,'\n');flush.console()
    command<-c('--vanilla','tools/benchmark-stable-scaling.R','worker',lib,fixtures,out,id,method)
    status<-system2(rs,vapply(command,shQuote,''))
    stopifnot(status==0L)
  }
  quit()
}
if(action=='worker') {
  id<-args[5];method<-args[6];stopifnot(method %in% c('maxLik','stable'))
  x<-load_case(id);affinity<-Sys.getenv('EMHAWKES_BENCH_AFFINITY','scheduler')
  if(affinity!='scheduler') {
    stopifnot(.Platform$OS.type=='windows',grepl('^[1-9][0-9]*$',affinity))
    actual<-system2('powershell.exe',c('-NoProfile','-NonInteractive','-Command',shQuote(sprintf(
      '$p = Get-Process -Id %d; $p.ProcessorAffinity = [IntPtr]%s; $p.ProcessorAffinity.ToInt64()',Sys.getpid(),affinity))),stdout=TRUE)
    stopifnot(length(actual)==1L,trimws(actual)==affinity)
  }
  warm<-x$arguments
  for(name in c('inter_arrival','type','mark'))warm[[name]]<-warm[[name]][seq_len(101L)]
  warm$control<-list(iterlim=2)
  invisible(suppressWarnings(suppressMessages(do.call(hfit,warm))))
  warnings<-messages<-character()
  invoke<-function()withCallingHandlers(
    do.call(hfit,c(x$arguments,list(final_hessian=method))),
    warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')},
    message=function(m){messages<<-c(messages,conditionMessage(m));invokeRestart('muffleMessage')})
  invisible(gc())
  elapsed<-system.time(fit<-invoke())
  diagnostics<-fit_diagnostics(fit)
  covariance<-tryCatch(suppressWarnings(vcov(fit)),error=function(e)
    matrix(NA_real_,length(coef(fit)),length(coef(fit)),dimnames=list(names(coef(fit)),names(coef(fit)))))
  se<-suppressWarnings(sqrt(diag(covariance)))
  public_summary<-suppressWarnings(summary(fit))$estimate
  stopifnot(isTRUE(all.equal(as.numeric(public_summary[,2]),as.numeric(se),tolerance=0)))
  stopifnot(diagnostics$computation$engine$used=='cpp',diagnostics$computation$gradient$used=='numeric')
  result<-list(id=id,method=method,meta=x$meta,value=fit[fields],covariance=covariance,
    se=se,diagnostics=diagnostics,information=information(fit$hessian),
    seconds=unname(elapsed[['elapsed']]),cpu_seconds=unname(elapsed[['user.self']]+elapsed[['sys.self']]),
    warnings=unique(warnings),messages=unique(messages),affinity=affinity)
  saveRDS(result,sample_path(id,method))
  cat('Result',id,method,'elapsed',result$seconds,
    'status',diagnostics$status,'nonfinite SE',sum(!is.finite(se)),'\n');flush.console()
  quit()
}
if(action=='compare') {
  rows<-coefficients<-steps<-notices<-list()
  for(id in design$id) {
    normal<-readRDS(sample_path(id,'maxLik'));stable<-readRDS(sample_path(id,'stable'))
    prior<-readRDS(file.path(fixtures,'samples',paste(id,'fit','numeric',1,'rds',sep='.')))
    retained<-setdiff(fields,'hessian')
    stopifnot(identical(normal$value,prior$value),
      identical(normal$value[retained],stable$value[retained]),
      identical(normal$value$hessian,stable$diagnostics$original_hessian),
      normal$value$code==0,stable$value$code==0)
    rows[[length(rows)+1L]]<-data.frame(id,events_used=normal$meta$events_used,
      dimensions=normal$meta$d,kernel_columns=normal$meta$k,parameters=length(normal$se),
      feedback=normal$meta$feedback,observed_baseline=normal$meta$observed,alpha_coupled=normal$meta$direct,
      seconds_default=normal$seconds,seconds_stable=stable$seconds,
      cpu_seconds_default=normal$cpu_seconds,cpu_seconds_stable=stable$cpu_seconds,
      additional_seconds=stable$seconds-normal$seconds,time_ratio=stable$seconds/normal$seconds,
      loglik_difference=abs(normal$value$maximum-stable$value$maximum),
      max_parameter_difference=max(abs(normal$value$estimate-stable$value$estimate)),
      code_default=normal$value$code,code_stable=stable$value$code,
      iterations_default=unname(normal$value$iterations),iterations_stable=unname(stable$value$iterations),
      nan_se_default=sum(is.nan(normal$se)),nan_se_stable=sum(is.nan(stable$se)),
      nonfinite_se_default=sum(!is.finite(normal$se)),nonfinite_se_stable=sum(!is.finite(stable$se)),
      min_scaled_information_default=normal$information$minimum,
      min_scaled_information_stable=stable$information$minimum,
      scaled_condition_stable=stable$information$condition,status_stable=stable$diagnostics$status,
      method_stable=stable$diagnostics$method,
      stable_step_levels=nrow(stable$diagnostics$checks),
      newton_decrement=stable$diagnostics$newton_decrement,
      maximum_standardized_score=stable$diagnostics$maximum_standardized_score)
    coefficients[[length(coefficients)+1L]]<-data.frame(id,parameter=names(normal$se),
      estimate_default=as.numeric(normal$value$estimate),estimate_stable=as.numeric(stable$value$estimate),
      se_default=as.numeric(normal$se),se_stable=as.numeric(stable$se))
    steps[[length(steps)+1L]]<-data.frame(id,stable$diagnostics$checks)
    for(z in list(normal,stable))for(kind in c('warnings','messages')) {
      for(note in z[[kind]])notices[[length(notices)+1L]]<-data.frame(id,method=z$method,kind,text=trimws(note))
    }
  }
  for(name in c('rows','coefficients','steps','notices')) {
    value<-get(name);value<-if(length(value))do.call(rbind,value) else data.frame()
    write.csv(value,file.path(out,paste0(if(name=='rows')'comparisons' else name,'.csv')),row.names=FALSE)
  }
  print(do.call(rbind,rows),row.names=FALSE);quit()
}
if(action=='reference') {
  rows<-coefficients<-list()
  for(id in design$id) {
    x<-load_case(id);a<-x$arguments;result<-readRDS(sample_path(id,'stable'))
    obs<-do.call(emhawkes:::normalize_observation,a[intersect(names(a),c('object','inter_arrival','type','mark'))])
    map<-emhawkes:::build_parameter_map(a$object,TRUE)
    objective<-emhawkes:::make_hawkes_objective(a$object,obs,a$lambda_component0,map,jacobian=x$jacobian)
    theta<-result$value$estimate;step<-2e-5
    H<-vapply(seq_along(theta),function(i) {
      plus<-minus<-theta;h<-step*max(1,abs(theta[i]))
      plus[i]<-plus[i]+h;minus[i]<-minus[i]-h
      (objective$gradient(plus)-objective$gradient(minus))/(2*h)
    },numeric(length(theta)))
    info<--(H+t(H))/2;C<-chol2inv(chol(info));se<-sqrt(diag(C))
    error<-max(abs(result$se-se)/pmax(abs(se),.Machine$double.eps))
    stopifnot(result$diagnostics$status=='ok',all(is.finite(result$se)),
      result$information$minimum>0,error<.02)
    normal<-readRDS(sample_path(id,'maxLik'))
    finite<-is.finite(normal$se)
    rows[[length(rows)+1L]]<-data.frame(id,step,stable_max_relative_se_error=error,
      stable_relative_covariance_error=norm(result$covariance-C,'F')/norm(C,'F'),
      default_finite_se_count=sum(finite),default_max_relative_error_among_finite=
        if(any(finite))max(abs(normal$se[finite]-se[finite])/se[finite]) else NA_real_)
    coefficients[[length(coefficients)+1L]]<-data.frame(id,parameter=names(theta),
      se_score_reference=se,se_stable=result$se,se_default=normal$se,
      stable_relative_error=abs(result$se-se)/se)
    cat('Verified score reference',id,'max relative SE error',error,'\n');flush.console()
    write.csv(do.call(rbind,rows),file.path(out,'reference-checks.csv'),row.names=FALSE)
    write.csv(do.call(rbind,coefficients),file.path(out,'reference-standard-errors.csv'),row.names=FALSE)
  }
  quit()
}
stop('Unknown action: ',action)
