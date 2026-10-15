# Rscript --vanilla tools/benchmark-gradient-scaling.R prepare|check|run|compare LIB OUT
# Individual workers: worker LIB OUT ID score|fit numeric|analytic ROUND
# ROUND=0 measures memory separately; positive rounds measure elapsed/CPU time.
args <- commandArgs(TRUE)
action <- args[1];lib <- args[2];out <- args[3]
.libPaths(c(lib,.libPaths()));library(emhawkes)
stopifnot(identical(normalizePath(dirname(find.package('emhawkes'))),normalizePath(lib)))
options(emhawkes.engine='auto')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source('tools/gradient-scaling-cases.R')
rs <- file.path(R.home('bin'),if(.Platform$OS.type=='windows')'Rscript.exe' else 'Rscript')
launch <- function(x) {status<-system2(rs,c('--vanilla',vapply(x,shQuote,'')));stopifnot(status==0L)}
if(action=='prepare') {
  grid <- data.frame(d=integer(),k=integer(),g=integer(),n=integer(),feedback=character(),observed=logical(),direct=logical(),fit=logical())
  add <- function(d,k,g=1L,n=10000L,feedback='diagonal',observed=TRUE,direct=FALSE,fit=FALSE) {
    grid[nrow(grid)+1L,]<<-list(d,k,g,n,feedback,observed,direct,fit)
  }
  for(n in c(1000L,10000L,100000L))add(2,4,n=n,fit=n!=100000L)
  for(d in c(1L,2L,4L,8L))add(d,d,fit=d!=2L)
  for(k in c(8L,16L))add(4,k,fit=k==16L)
  for(g in c(1L,2L,4L,8L))add(8,16,g,fit=g %in% c(1L,8L))
  for(g in c(2L,4L))add(4,8,g,fit=g==4L)
  add(4,8,4,feedback='dense')
  add(8,16,8,feedback='dense',fit=TRUE)
  add(8,16,8,observed=FALSE,fit=TRUE)
  add(4,8,4,feedback='dense',direct=TRUE,fit=TRUE)
  grid$id <- with(grid,sprintf('n%d-d%d-k%d-p%d-%s-%s%s',n,d,k,3*g,feedback,
    ifelse(observed,'observed','constant'),ifelse(direct,'-coupled','')))
  stopifnot(!anyDuplicated(grid$id))
  dir.create(file.path(out,'fixtures'),showWarnings=FALSE)
  # Each parameter-sharing variant uses the same generating model and data.
  family <- with(grid,paste(d,k,feedback,observed,direct,sep=':'))
  for(key in unique(family)) {
    indices <- which(family==key);row <- grid[indices[1],]
    truth <- scaling_markov_model(row$d,row$k,1L,row$feedback,row$observed,row$direct)
    seed <- 75100L+row$d*100L+row$k+1000L*(row$feedback=='dense')+3000L*(!row$observed)+5000L*row$direct
    set.seed(seed)
    path <- simulate_scaling_markov(truth,min(10000L,max(grid$n[indices])))
    for(i in indices) {
      meta <- grid[i,];model <- scaling_markov_model(meta$d,meta$k,meta$g,meta$feedback,meta$observed,meta$direct)
      map <- emhawkes:::build_parameter_map(model$object,TRUE)
      p <- map$start+c(rep(.08,meta$g),rep(-.1,meta$g),rep(.07,meta$g))
      v <- emhawkes:::restore_parameters(map,p)
      for(name in names(v))slot(model$object,name)<-v[[name]]
      use <- seq_len(meta$n+1L)
      # The 100k case tests evaluation cost on a valid observed sequence,
      # not parameter recovery. Full-fit cases use the simulated paths above.
      data <- if(meta$n>10000L) {
        set.seed(seed+99000L)
        list(inter_arrival=c(0,rexp(meta$n,rate=meta$d*.7)),
          type=c(0L,sample.int(meta$d,meta$n,replace=TRUE)),mark=c(0L,rgeom(meta$n,.7)+1L))
      } else path
      fixture <- list(meta=meta,jacobian=model$jacobian,
        arguments=list(object=model$object,inter_arrival=data$inter_arrival[use],type=data$type[use],
          mark=data$mark[use],lambda_component0=model$initial,method='BFGS',control=list(iterlim=500,reltol=1e-10)))
      saveRDS(fixture,file.path(out,'fixtures',paste0(meta$id,'.rds')))
    }
    cat('Prepared',key,'events',max(grid$n[indices]),'\n');flush.console()
  }
  write.csv(grid,file.path(out,'design.csv'),row.names=FALSE)
  writeLines(c(paste('source_commit',system2('git',c('rev-parse','HEAD'),stdout=TRUE)),
    paste('R',getRversion()),paste('maxLik',packageVersion('maxLik')),paste('Rcpp',packageVersion('Rcpp'))),file.path(out,'manifest.txt'))
  writeLines(trimws(capture.output(sessionInfo()),which='right'),file.path(out,'session-info.txt'))
  quit()
}
grid <- read.csv(file.path(out,'design.csv'),stringsAsFactors=FALSE)
# Repeat the near-break-even full fit to distinguish timing noise from a regression.
fit_rounds <- function(id)if(id=='n10000-d8-k16-p3-diagonal-observed')1:3 else 1L
if(action=='run') {
  dir.create(file.path(out,'samples'),showWarnings=FALSE)
  for(task in c('score','fit'))for(id in grid$id[if(task=='fit')grid$fit else rep(TRUE,nrow(grid))])
    for(round in if(task=='fit')fit_rounds(id) else c(1L,2L,3L,0L))for(method in if(round %% 2)c('numeric','analytic') else c('analytic','numeric')) {
      file <- file.path(out,'samples',paste(id,task,method,round,'rds',sep='.'))
      if(file.exists(file))next
      cat(id,task,method,round,'\n');flush.console()
      launch(c('tools/benchmark-gradient-scaling.R','worker',lib,out,id,task,method,round))
    }
  quit()
}
fixture <- function(id)readRDS(file.path(out,'fixtures',paste0(id,'.rds')))
setup <- function(x,method) {
  a<-x$arguments
  obs<-do.call(emhawkes:::normalize_observation,a[intersect(names(a),c('object','inter_arrival','type','mark'))])
  map<-emhawkes:::build_parameter_map(a$object,TRUE)
  obj<-emhawkes:::make_hawkes_objective(a$object,obs,a$lambda_component0,map,
    derivatives=method=='analytic',jacobian=if(method=='analytic')x$jacobian else NULL)
  list(obs=obs,map=map,objective=obj)
}
central <- function(f,p,step) vapply(seq_along(p),function(i) {
  a<-b<-p;h<-step*max(1,abs(p[i]));a[i]<-a[i]+h;b[i]<-b[i]-h;(f(a)-f(b))/(2*h)
},0)
if(action=='check') {
  rows<-list()
  for(id in grid$id) {
    x<-fixture(id);a<-setup(x,'analytic');b<-setup(x,'numeric');p<-a$map$start
    reference<-function(q) {
      v<-emhawkes:::restore_parameters(a$map,q);h<-x$arguments$object
      for(name in names(v))slot(h,name)<-v[[name]]
      old<-options(emhawkes.engine='r');on.exit(options(old))
      as.numeric(suppressMessages(logLik(h,x$arguments$inter_arrival,x$arguments$type,x$arguments$mark,
        lambda_component0=x$arguments$lambda_component0)))
    }
    direction<-cos(seq_along(p));direction<-direction/sqrt(sum(direction^2))
    for(point in c('start','nearby')) {
      q<-if(point=='start')p else p+.03*direction
      ll<-as.numeric(a$objective$loglik(q));r<-reference(q)
      stopifnot(identical(a$objective$loglik(q),b$objective$loglik(q)),abs(ll-r)<1e-7)
      g<-a$objective$gradient(q)
      for(step in c(1e-4,2e-5)) {
        fd<-central(b$objective$loglik,q,step)
        error<-max(abs(g-fd)/pmax(1,abs(fd)))
        stopifnot(error<5e-5)
        rows[[length(rows)+1L]]<-cbind(x$meta,data.frame(point,step,loglik_R=r,loglik_cpp=ll,
          loglik_abs_difference=abs(ll-r),gradient_max_abs_difference=max(abs(g-fd)),gradient_max_scaled_difference=error))
      }
      dir_fd<-(reference(q+2e-5*direction)-reference(q-2e-5*direction))/4e-5
      stopifnot(abs(sum(g*direction)-dir_fd)/max(1,abs(dir_fd))<5e-5)
    }
    cat('Checked',id,'\n');flush.console()
    write.csv(do.call(rbind,rows),file.path(out,'derivative-checks.csv'),row.names=FALSE)
  }
  quit()
}
if(action=='worker') {
  id<-args[4];task<-args[5];method<-args[6];round<-as.integer(args[7]);x<-fixture(id)
  x$meta$events_used <- if(task=='fit' && (x$meta$d>=4L || x$meta$g>1L))min(2500L,x$meta$n) else x$meta$n
  if(task=='fit')for(name in c('inter_arrival','type','mark'))
    x$arguments[[name]]<-x$arguments[[name]][seq_len(x$meta$events_used+1L)]
  affinity<-Sys.getenv('EMHAWKES_BENCH_AFFINITY','scheduler')
  if(affinity!='scheduler') {
    stopifnot(.Platform$OS.type=='windows',grepl('^[1-9][0-9]*$',affinity))
    actual<-system2('powershell.exe',c('-NoProfile','-NonInteractive','-Command',shQuote(sprintf(
      '$p = Get-Process -Id %d; $p.ProcessorAffinity = [IntPtr]%s; $p.ProcessorAffinity.ToInt64()',Sys.getpid(),affinity))),stdout=TRUE)
    stopifnot(length(actual)==1L,trimws(actual)==affinity)
  }
  # Warm both callback/optimizer paths on a small prefix, never the measured cache.
  warm<-x
  for(name in c('inter_arrival','type','mark'))warm$arguments[[name]]<-warm$arguments[[name]][seq_len(min(101L,x$meta$n+1L))]
  if(task=='fit') {
    warm$arguments$control<-list(iterlim=2)
    invisible(suppressWarnings(suppressMessages(do.call(hfit,c(warm$arguments,
      if(method=='analytic')list(gradient='analytic',jacobian=x$jacobian) else list())))))
  } else {
    w<-setup(warm,method);q<-w$map$start
    invisible(w$objective$loglik(q))
    invisible(if(method=='analytic')w$objective$gradient(q) else maxLik::numericGradient(w$objective$loglik,q))
  }
  if(task=='score') {s<-setup(x,method);p<-s$map$start;counter<-0L}
  warnings<-character()
  operation<-if(task=='fit')function()suppressMessages(do.call(hfit,c(x$arguments,
    if(method=='analytic')list(gradient='analytic',jacobian=x$jacobian) else list()))) else function() {
    counter<<-counter+1L;q<-p;q[1]<-q[1]+1e-7*counter
    ll<-s$objective$loglik(q)
    g<-if(method=='analytic')s$objective$gradient(q) else
      setNames(as.numeric(maxLik::numericGradient(s$objective$loglik,q)),names(p))
    list(theta=q,loglik=ll,gradient=g)
  }
  invoke<-function()withCallingHandlers(operation(),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  rss<-function()if(.Platform$OS.type=='windows')as.numeric(system2('powershell.exe',c('-NoProfile','-NonInteractive','-Command',
    shQuote(sprintf('(Get-Process -Id %d).PeakWorkingSet64',Sys.getpid()))),stdout=TRUE)) else NA_real_
  gc();memory<-NULL
  if(task=='fit') {fit_before<-gc(reset=TRUE);fit_rss_before<-rss()}
  if(round==0L) {
    before<-gc(reset=TRUE);rss_before<-rss();profile<-tempfile()
    # Full-fit heap/RSS peaks are measured without a potentially enormous trace.
    if(task=='score')Rprofmem(profile)
    value<-tryCatch(invoke(),finally=if(task=='score')Rprofmem(NULL))
    after<-gc();rss_after<-rss();bytes<-allocations<-largest<-0
    if(task=='score') {
      con<-file(profile,'r')
      repeat {lines<-readLines(con,n=10000L,warn=FALSE);if(!length(lines))break
        z<-suppressWarnings(as.numeric(sub(' .*','',substr(lines,1,24))));z<-z[is.finite(z)]
        bytes<-bytes+sum(z);allocations<-allocations+length(z);if(length(z))largest<-max(largest,z)}
      close(con);unlink(profile)
    }
    memory<-c(allocated_bytes=if(task=='score')bytes else NA_real_,allocation_count=if(task=='score')allocations else NA_real_,
      largest_block_bytes=if(task=='score')largest else NA_real_,baseline_r_heap_mib=sum(before[,2]),
      peak_r_heap_mib=sum(after[,6]),incremental_r_heap_mib=sum(after[,6])-sum(before[,2]),
      baseline_peak_rss_bytes=rss_before,peak_rss_bytes=rss_after)
    batch<-1L;seconds<-cpu<-NA_real_
  } else {
    batch<-if(task=='fit')1L else if(x$meta$n<=1000)10L else 1L
    t<-system.time(for(i in seq_len(batch))value<-invoke())
    seconds<-unname(t[['elapsed']])/batch;cpu<-unname(t[['user.self']]+t[['sys.self']])/batch
  }
  if(task=='fit') {
    fit_after<-gc();fit_rss_after<-rss()
    memory<-c(allocated_bytes=NA_real_,allocation_count=NA_real_,largest_block_bytes=NA_real_,
      baseline_r_heap_mib=sum(fit_before[,2]),peak_r_heap_mib=sum(fit_after[,6]),
      incremental_r_heap_mib=sum(fit_after[,6])-sum(fit_before[,2]),
      baseline_peak_rss_bytes=fit_rss_before,peak_rss_bytes=fit_rss_after)
    computation<-fit_diagnostics(value)$computation
    stopifnot(identical(computation$engine$used,'cpp'),computation$gradient$used==method)
    value<-value[c('estimate','maximum','gradient','hessian','code','iterations')]
  } else computation<-NULL
  saveRDS(list(meta=x$meta,id=id,task=task,method=method,round=round,batch=batch,
    seconds=seconds,cpu_seconds=cpu,memory=memory,value=value,computation=computation,
    warnings=unique(warnings),affinity=affinity),file.path(out,'samples',paste(id,task,method,round,'rds',sep='.')))
  cat('Result',id,task,method,round,'seconds',seconds,'\n');flush.console()
  quit()
}
if(action=='compare') {
  samples<-lapply(list.files(file.path(out,'samples'),pattern='[.]rds$',full.names=TRUE),readRDS)
  timings<-memory<-checks<-estimates<-methods<-list()
  for(task in c('score','fit'))for(id in grid$id[if(task=='fit')grid$fit else rep(TRUE,nrow(grid))]) {
    selected<-Filter(function(x)x$id==id && x$task==task,samples)
    stopifnot(length(selected)==if(task=='fit')2L*length(fit_rounds(id)) else 8L)
    ref<-Filter(function(x)x$method=='numeric' && x$round==1L,selected)[[1]]
    for(x in selected) {
      stopifnot(length(x$warnings)==0L)
      if(task=='score') {
        # Round zero has one proposal, so its numeric reference uses that proposal too.
        other<-Filter(function(z)z$method=='numeric' && z$round==x$round,selected)[[1]]
        stopifnot(identical(x$value$theta,other$value$theta),identical(x$value$loglik,other$value$loglik))
        err<-max(abs(x$value$gradient-other$value$gradient)/pmax(1,abs(other$value$gradient)))
        stopifnot(err<2e-4)
        checks[[length(checks)+1L]]<-cbind(x$meta,data.frame(task,method=x$method,round=x$round,
          loglik_difference=0,max_parameter_difference=NA_real_,max_scaled_parameter_difference=NA_real_,
          max_scaled_gradient_difference=err,code=NA_integer_,iterations=NA_integer_))
      } else {
        repeat_ref<-Filter(function(z)z$method==x$method && z$round==1L,selected)[[1]]
        stopifnot(isTRUE(all.equal(x$value,repeat_ref$value,tolerance=0)))
        computation<-x$computation
        stopifnot(computation$hessian$used=='maxLik_numeric')
        methods[[length(methods)+1L]]<-data.frame(id,method=x$method,round=x$round,
          engine_requested=computation$engine$requested,engine_used=paste(computation$engine$used,collapse='+'),
          gradient_requested=computation$gradient$requested,gradient_used=computation$gradient$used,
          hessian_requested=computation$hessian$requested,hessian_used=computation$hessian$used)
        ll<-abs(x$value$maximum-ref$value$maximum)
        pd<-max(abs(x$value$estimate-ref$value$estimate))
        spd<-max(abs(x$value$estimate-ref$value$estimate)/pmax(1,abs(ref$value$estimate)))
        stopifnot(x$value$code==0,ll<1e-5,spd<2e-3)
        checks[[length(checks)+1L]]<-cbind(x$meta,data.frame(task,method=x$method,round=x$round,
          loglik_difference=ll,max_parameter_difference=pd,max_scaled_parameter_difference=spd,
          max_scaled_gradient_difference=NA_real_,code=x$value$code,iterations=unname(x$value$iterations)))
        if(x$round==1L)estimates[[length(estimates)+1L]]<-data.frame(id,method=x$method,
          parameter=names(x$value$estimate),estimate=as.numeric(x$value$estimate),gradient=as.numeric(x$value$gradient))
      }
      if(x$round>0L)timings[[length(timings)+1L]]<-cbind(x$meta,data.frame(task,method=x$method,round=x$round,
        batch=x$batch,seconds=x$seconds,cpu_seconds=x$cpu_seconds,affinity=x$affinity))
      if(!is.null(x$memory))memory[[length(memory)+1L]]<-cbind(x$meta,data.frame(task,method=x$method,round=x$round,as.list(x$memory)))
    }
  }
  timings<-do.call(rbind,timings)
  summary<-aggregate(seconds~id+task+method,timings,median)
  pairs<-merge(subset(summary,method=='numeric',select=-method),subset(summary,method=='analytic',select=-method),
    by=c('id','task'),suffixes=c('_numeric','_analytic'))
  pairs$speedup<-pairs$seconds_numeric/pairs$seconds_analytic
  for(name in c('timings','summary','pairs','memory','checks','estimates','methods')) {
    value<-get(name);if(is.list(value)&&!is.data.frame(value))value<-do.call(rbind,value)
    write.csv(value,file.path(out,paste0(name,'.csv')),row.names=FALSE)
  }
  print(pairs,row.names=FALSE);quit()
}
stop('Unknown action: ',action)
