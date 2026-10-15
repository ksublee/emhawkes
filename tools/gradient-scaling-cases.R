# Markov fixtures for derivative scaling. All state access is current-component.
scaling_markov_model <- function(dimensions, kernels, groups = 1L,
    feedback = c('diagonal', 'dense'), observed_mu = TRUE, direct_alpha = FALSE) {
  feedback <- match.arg(feedback)
  d <- as.integer(dimensions); k <- as.integer(kernels); g <- as.integer(groups)
  stopifnot(d >= 1L, k %% d == 0L, d %% g == 0L)
  m <- d * k
  group <- rep(seq_len(g), each=d/g)
  entry_group <- rep(group,k)
  mu_names <- paste0('log_mu',seq_len(g))
  alpha_names <- paste0('log_alpha',seq_len(g))
  beta_names <- paste0('log_beta',seq_len(g))
  starts <- list(mu=setNames(rep(log(.42),g),mu_names),
    alpha=setNames(rep(log(.55),g),alpha_names),
    beta=setNames(rep(0,g),beta_names))
  weight <- matrix(0,d,k)
  for(i in seq_len(d))for(j in seq_len(k)) {
    same <- ((j-1L) %% d)+1L == i
    weight[i,j] <- if(d==1L)1/k else if(same).7/(k/d) else .3/((d-1)*(k/d))
  }
  rates <- outer(seq(.95,1.05,length.out=d),seq(.9,2.1,length.out=k))
  mu_weight <- seq(.9,1.1,length.out=d)
  q <- .12
  counts_factor <- function(n,N) 1+.08*(N[max(1L,n-1L),] %% 3)
  mu <- if(observed_mu) function(param,n,N,...) {
    matrix(mu_weight*exp(param[group])*counts_factor(n,N),d)
  } else function(param) matrix(mu_weight*exp(param[group]),d)
  formals(mu)$param <- starts$mu
  alpha <- function(param) weight*matrix(exp(param[entry_group]),d,k)
  beta <- function(param) rates*matrix(exp(param[entry_group]),d,k)
  formals(alpha)$param <- starts$alpha;formals(beta)$param <- starts$beta
  impact <- function(n,mark,alpha,lambda_component_n,...) {
    factor <- 1+.1*(mark[n]-1)
    c <- lambda_component_n
    extra <- if(feedback=='diagonal') q*c/(1+c) else matrix(q*mean(c)/(1+mean(c)),d,k)
    if(direct_alpha)extra <- extra+.1*alpha
    factor*extra
  }
  local_jac <- function(param,values) {
    J <- matrix(0,length(values),length(param),dimnames=list(NULL,names(param)))
    index <- if(length(values)==d)group else entry_group
    J[cbind(seq_along(values),index)] <- as.vector(values)
    J
  }
  mu_jac <- if(observed_mu) function(param,n,N,...) {
    J <- matrix(0,d,length(param),dimnames=list(NULL,names(param)))
    J[cbind(seq_len(d),match(mu_names[group],names(param)))] <-
      mu_weight*exp(param[mu_names[group]])*counts_factor(n,N)
    J
  } else function(param) local_jac(param,matrix(mu_weight*exp(param[group]),d))
  alpha_jac <- function(param) local_jac(param,alpha(param))
  beta_jac <- function(param) local_jac(param,beta(param))
  impact_param <- function(param,n,mark,alpha,...) {
    J <- matrix(0,m,length(param),dimnames=list(NULL,names(param)))
    if(direct_alpha) J[cbind(seq_len(m),match(alpha_names[entry_group],names(param)))] <-
      .1*(1+.1*(mark[n]-1))*as.vector(alpha)
    J
  }
  impact_state <- function(n,mark,lambda_component_n,...) {
    factor <- q*(1+.1*(mark[n]-1))
    c <- as.vector(lambda_component_n)
    if(feedback=='diagonal') diag(factor/(1+c)^2,m) else
      matrix(factor/(m*(1+mean(c))^2),m,m)
  }
  h <- new('hspec',mu=mu,alpha=alpha,beta=beta,impact=impact,
    rmark=function()rgeom(1,.7)+1,
    type_col_map=lapply(seq_len(d),function(i)seq(i,k,by=d)))
  list(object=h,jacobian=list(mu=mu_jac,alpha=alpha_jac,beta=beta_jac,
      impact=list(param=if(direct_alpha)impact_param else
        matrix(0,m,3L*g,dimnames=list(NULL,c(mu_names,alpha_names,beta_names))),state=impact_state)),
    initial=matrix(.03/k,d,k),dimensions=d,kernels=k,parameters=3L*g,
    groups=g,feedback=feedback,observed_mu=observed_mu,direct_alpha=direct_alpha)
}

# These fixtures use only current components and cumulative counts. Simulate
# bounded chunks with the existing hsim(), carrying the post-event Markov state.
# This avoids including a large simulation history in benchmark preparation.
simulate_scaling_markov <- function(model,events,chunk_size=500L) {
  dt<-mark<-numeric(events+1L);type<-integer(events+1L)
  state<-model$initial;N0<-Nc0<-rep(0,model$dimensions);done<-0L;checkpoints<-list()
  while(done<events) {
    count<-min(chunk_size,events-done)
    path<-hsim(model$object,size=count+1L,lambda_component0=state,N0=N0,Nc0=Nc0)
    stopifnot(length(path$inter_arrival)==count+1L)
    use<-seq_len(count)+1L;target<-done+use
    dt[target]<-path$inter_arrival[use];type[target]<-path$type[use];mark[target]<-path$mark[use]
    state<-matrix(path$rambda_component[count+1L,],model$dimensions,model$kernels,byrow=TRUE)
    N0<-as.numeric(path$N[count+1L,]);Nc0<-as.numeric(path$Nc[count+1L,])
    done<-done+count
    checkpoints[[length(checkpoints)+1L]]<-list(row=done+1L,state=as.vector(t(state)))
  }
  # Replay the full joined path once to verify that no state was lost at joins.
  replay<-suppressMessages(logLik(model$object,dt,type,mark,lambda_component0=model$initial,infer=TRUE))
  for(x in checkpoints)stopifnot(isTRUE(all.equal(as.numeric(replay$rambda_component[x$row,]),x$state,tolerance=1e-11)))
  list(inter_arrival=dt,type=type,mark=mark)
}
