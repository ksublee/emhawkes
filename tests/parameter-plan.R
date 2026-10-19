library(emhawkes)
fails <- function(expr,text) {
 error <- tryCatch(force(expr),error=identity)
 stopifnot(inherits(error,"error"),grepl(text,conditionMessage(error),fixed=TRUE))
}
set.seed(523);saved <- .Random.seed
# Explicit vectors retain distinct names even when their values are equal.
plan <- parameter_plan(list(mu=matrix(c(.2,.2),2),
 alpha=matrix(c(.1,.3,.1,.4),2,byrow=TRUE),residual=c(shape=2,rate=2)),named="residual")
stopifnot(identical(plan$start,c(mu1=.2,alpha1.1=.1,alpha1.2=.3,alpha2.2=.4,shape=2,rate=2)),
 environmentIsLocked(plan),bindingIsLocked("start",plan),bindingIsLocked("restore",plan))
empty <- parameter_plan(list(optional=NULL))
stopifnot(identical(empty$start,numeric()),identical(empty$restore(empty$start),list(optional=NULL)))
values <- plan$restore(plan$start)
stopifnot(identical(values$mu,matrix(.2,2)),
 identical(values$alpha,matrix(c(.1,.3,.1,.4),2,byrow=TRUE)),
 identical(values$residual,c(shape=2,rate=2)),
 identical(values,plan$restore(rev(plan$start))),
 identical(values,plan$restore(unname(plan$start))))
full <- parameter_plan(list(alpha=matrix(c(.1,.3,.1,.4),2,byrow=TRUE)),reduced=FALSE)
stopifnot(identical(full$start,c(alpha1.1=.1,alpha2.1=.1,alpha1.2=.3,alpha2.2=.4)))
# Shared coordinates cross component boundaries; preserve first-value warning.
warnings <- character()
shared <- withCallingHandlers(parameter_plan(list(mu=matrix(.2),
 residual=c(mu1=.7,scale=1)),named="residual"),warning=function(w) {
 warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
})
stopifnot(length(warnings)==1L,grepl("using the first value for: mu1",warnings,fixed=TRUE),
 identical(shared$start,c(mu1=.2,scale=1)),
 identical(shared$restore(c(mu1=.4,scale=2))$residual,c(mu1=.4,scale=2)))
# Numeric inputs are captured; function environments/default binding stay intact.
original <- matrix(.2)
fixed <- parameter_plan(list(mu=original));original[] <- 9
stopifnot(identical(fixed$restore(fixed$start)$mu,matrix(.2)))
local({
 multiplier <- 2; calls <- 0L
 callback <- function(param=c(a=.1), n=1) {calls <<- calls+1L;matrix(multiplier*param["a"]*n)}
 frozen <- formals(callback)
 p <- parameter_plan(list(alpha=callback,optional=NULL))
 stopifnot(calls==0L,identical(formals(callback),frozen))
 bound <- p$restore(c(a=.4))$alpha
 stopifnot(identical(bound(n=2),matrix(1.6)),calls==1L,
   identical(formals(callback),frozen))
 multiplier <- 3
 stopifnot(identical(bound(n=2),matrix(3*.4*2)))
})
restored <- unserialize(serialize(plan,NULL))
stopifnot(identical(restored$start,plan$start),
 identical(restored$restore(plan$start),values),identical(saved,.Random.seed))
# The plan preserves numerical proposals; the model decides their admissibility.
proposal <- plan$start;proposal[1] <- Inf
stopifnot(all(is.infinite(plan$restore(proposal)$mu)))
fails(parameter_plan(list(.1)),"components")
fails(parameter_plan(list(x=1),named="missing"),"named")
fails(parameter_plan(list(x=c(a=1,a=2)),named="x"),"unique")
fails(parameter_plan(list(x=c(1,2))),"matrix")
fails(parameter_plan(list(x=1),reduced=NA),"reduced")
fails(plan$restore(c(mu1=.2)),"full parameter vector")
bad <- plan$start;names(bad)[1] <- "other"
fails(plan$restore(bad),"exactly its names")
fails(plan$restore(matrix(plan$start,1)),"full parameter vector")
fails({plan$start <- numeric()},"locked")
stopifnot(identical(saved,.Random.seed))
cat("Parameter coordinates, sharing, vector restoration, callbacks, capture, serialization and failure contracts passed.\n")
