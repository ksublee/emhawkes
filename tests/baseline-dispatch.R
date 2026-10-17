# Legacy optimizer checks require the optional maxLik package.
if (requireNamespace("maxLik", quietly = TRUE)) {
library(emhawkes)
old <- options(emhawkes.engine="auto")
jit <- compiler::enableJIT(0)
collect <- emhawkes:::hawkes_baseline_cpp
fails <- function(expr, pattern) {
 error <- tryCatch(force(expr),error=identity)
 stopifnot(inherits(error,"error"),grepl(pattern,conditionMessage(error)))
 invisible(error)
}
# Saved lazy arguments retain every event index and the original shared inputs.
saved <- list()
callback <- function(n, Nc, alpha, param=c(rate=.3), extra=stop("unused default")) {
 default <- formals(sys.function())$param
 saved[[length(saved)+1L]] <<- function() list(n=n,Nc=Nc,alpha=alpha,param=param,default=default)
 value <- matrix(param + .01 * Nc[max(1,n-1),1])
 Nc[1,1] <- 999; alpha[1] <- 999
 invisible(gc(FALSE))
 value
}
arguments <- list(n=1L,Nc=matrix(0:4),alpha=matrix(.1))
before <- serialize(arguments,NULL)
for(integer_n in c(FALSE,TRUE)) for(infer in c(FALSE,TRUE)) for(initial in list(NULL,matrix(.3))) {
 saved <- list()
 value <- collect(callback,arguments,1L,5L,initial,infer,integer_n)
 snapshots <- lapply(saved,function(f)f())
 indices <- if(is.null(initial))1:6 else 2:6
 stopifnot(identical(vapply(snapshots,`[[`,0,"n"),as.numeric(indices)))
 stopifnot(all(vapply(snapshots,function(x)typeof(x$n)==if(integer_n)"integer" else "double",TRUE)))
 saved <- list()
 for(n in indices) { a <- arguments; a$n <- if(integer_n)as.integer(n) else as.double(n); do.call(callback,a) }
 expected_snapshots <- lapply(saved,function(f)f())
 stopifnot(identical(snapshots,expected_snapshots))
 expected <- .3 + .01 * c(0,0:4)
 stopifnot(identical(value$path,matrix(expected[seq_len(if(infer)6L else 5L)],nrow=1)),
   identical(before,serialize(arguments,NULL)))
}
# Functions with no runtime formals are valid observed baselines too.
stopifnot(identical(collect(function()matrix(.3),list(),1L,3L)$path,matrix(.3,1,3)))
for(fun in list(function()matrix(.3),function(param=c(rate=.3),extra=0)matrix(param+extra))) {
 h0 <- new("hspec",mu=fun,alpha=.1,beta=1.2)
 options(emhawkes.engine="r");ref0<-logLik(h0,c(0,.2,.3),lambda_component0=.1,infer=TRUE)
 options(emhawkes.engine="cpp");stopifnot(identical(logLik(h0,c(0,.2,.3),lambda_component0=.1,infer=TRUE),ref0))
 m0<-emhawkes:::build_parameter_map(h0)
 obj0<-emhawkes:::make_hawkes_objective(h0,emhawkes:::normalize_observation(h0,c(0,.2,.3)),.1,m0,derivatives=FALSE)
 stopifnot(isTRUE(all.equal(as.numeric(obj0$loglik(m0$start)),as.numeric(ref0$loglikelihood),tolerance=1e-12)))
}
# Native calls are reentrant and preserve user condition classes and defaults.
inner <- function(n)matrix(.2)
outer <- function(n) {
 stopifnot(identical(collect(inner,list(n=1L),1L,2L)$path,matrix(.2,1,2)))
 matrix(.3)
}
stopifnot(identical(collect(outer,list(n=1L),1L,3L)$path,matrix(.3,1,3)))
custom <- function(n) { if(n==3)stop(structure(list(message="baseline sentinel",call=NULL),class=c("baseline_sentinel","error","condition")));matrix(.3) }
e <- fails(collect(custom,list(n=1L),1L,4L),"baseline sentinel")
stopifnot(inherits(e,"baseline_sentinel"))
for(value in list(1,matrix(.3,2),structure(matrix(.3),class=c("custom","matrix","array")))) {
 fun <- local({v<-value;function(n)v})
 stopifnot(!is.null(collect(fun,list(n=1L),1L,3L)$unsupported))
}
for(value in c(-1,Inf,NaN)) {
 fun <- local({v<-value;function(n)matrix(v)})
 stopifnot(isTRUE(collect(fun,list(n=1L),1L,3L)$invalid))
}
# Only the terminal value is unused without inference; terminal errors propagate.
terminal <- function(n)if(n==4L)NULL else matrix(1L)
stopifnot(identical(collect(terminal,list(n=1L),1L,3L)$path,matrix(1,1,3)),
 !is.null(collect(terminal,list(n=1L),1L,3L,infer=TRUE)$unsupported))
fails(collect(inner,list(n=1L),1L,0L),"dimensions")
fails(collect(inner,list(1L),1L,2L),"named")
fails(collect(inner,list(n=1L),1L,2L,initial=matrix(0,2)),"initial baseline")
# A latent baseline read discovered after initialization keeps its typed fallback.
h <- new("hspec",mu=function(param=c(rate=.3),n,...) {
 if(n>=3)matrix(param + 0*sum(list(...)$lambda_component_n)) else matrix(param)
},alpha=.1,beta=1.2)
dt <- c(0,.2,.3,.4)
options(emhawkes.engine="r"); ref <- logLik(h,dt,lambda_component0=.1,infer=TRUE)
options(emhawkes.engine="auto"); candidate <- logLik(h,dt,lambda_component0=.1,infer=TRUE)
stopifnot(identical(candidate,ref))
options(emhawkes.engine="cpp")
fails(logLik(h,dt,lambda_component0=.1),"latent intensit")
# A fit also falls back once when a capability failure occurs inside the batch.
fallbacks <- 0L
unsupported <- emhawkes:::stop_cpp_unsupported
assignInNamespace("stop_cpp_unsupported",function(message){fallbacks<<-fallbacks+1L;unsupported(message)},"emhawkes")
options(emhawkes.engine="auto")
auto <- suppressWarnings(suppressMessages(hfit(h,dt,lambda_component0=.1,control=list(iterlim=2))))
assignInNamespace("stop_cpp_unsupported",unsupported,"emhawkes")
stopifnot(fallbacks==1L)
options(emhawkes.engine="r")
r <- suppressWarnings(suppressMessages(hfit(h,dt,lambda_component0=.1,control=list(iterlim=2))))
for(field in c("estimate","maximum","gradient","hessian","code","iterations"))stopifnot(identical(auto[[field]],r[[field]]))
options(old);invisible(compiler::enableJIT(jit))

}
