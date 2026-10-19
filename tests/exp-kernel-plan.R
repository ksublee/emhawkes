library(emhawkes)
fails <- function(expr, text) {
  error <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(error, "error"), grepl(text, conditionMessage(error), fixed = TRUE))
}
set.seed(991)
saved <- .Random.seed
comparisons <- 0L
for (d in 1:3) for (multiple in 1:2) for (marked in c(FALSE, TRUE)) {
  k <- d * multiple
  kernel <- exp_kernel_spec(matrix(seq_len(d*k)/100, d),
    matrix(rep(seq_len(d)/2+1, k), d),
    eta = if (marked) matrix(.02, d, k) else NULL,
    type_col_map = lapply(seq_len(d), function(i) rev(unique(c(i,k)))))
  mu <- seq_len(d)/10
  plan <- exp_kernel_plan(kernel, mu)
  stopifnot(environmentIsLocked(plan), bindingIsLocked("step", plan),
    bindingIsLocked("wait", plan), identical(sort(ls(plan)), c("step","wait")))
  component <- matrix(seq_len(d*k)/20, d)
  component[1,1] <- -.02
  for (dt in c(0, 1e-14, .2, 3)) for (type in c(list(NULL),as.list(seq_len(d))))
    for (mark in c(.5,1,2)) for (integrate in c(FALSE,TRUE)) {
      reference <- exp_kernel_step(kernel,component,dt,type,mark,integrate)
      stopifnot(identical(reference,plan$step(component,dt,type,mark,integrate)))
      comparisons <- comparisons + 1L
    }
  if (k == 1L) component <- abs(component)
  for (target in list(rep(0,d),rep(.4,d),rep(Inf,d))) {
    reference <- exp_kernel_wait(kernel,component,mu,target)
    stopifnot(identical(reference,plan$wait(component,target)),
      identical(reference,plan$wait(component,target)))
  }
  frozen_kernel <- kernel; frozen_mu <- mu
  kernel$alpha[] <- 99; kernel$beta[] <- 99; kernel$type_col_map <- NULL
  mu[] <- 99
  stopifnot(identical(plan$step(component,.2,1,2),
    exp_kernel_step(frozen_kernel,component,.2,1,2)),
    identical(plan$wait(component,rep(.4,d)),
    exp_kernel_wait(frozen_kernel,component,frozen_mu,rep(.4,d))))
  restored <- unserialize(serialize(plan,NULL))
  stopifnot(identical(plan$step(component,.2,1),restored$step(component,.2,1)),
    identical(plan$wait(component,rep(.4,d)),restored$wait(component,rep(.4,d))))
}
kernel <- exp_kernel_spec(matrix(.1,2,2),matrix(c(1,2,1,2),2))
plan <- exp_kernel_plan(kernel,0)
stopifnot(identical(plan$wait(0,c(0,1)),c(0,Inf)),
  identical(plan$wait(.1,c(Inf,.1)),c(Inf,Inf)))
fails(exp_kernel_plan(kernel,-1),"mu")
fails(exp_kernel_plan(kernel,NA_real_),"mu")
fails(plan$step(numeric(0),0),"component")
fails(plan$step(matrix(0,1,4),0),"component")
fails(plan$step(Inf,0),"component")
fails(plan$step(0,-1),"dt")
fails(plan$step(0,0,3),"type")
fails(plan$step(0,0,1,NA_real_),"mark")
fails(plan$step(0,0,integrate=NA),"integrate")
fails(plan$wait(matrix(0,1,4),c(1,1)),"component")
fails(plan$wait(NA_real_,c(1,1)),"component")
fails(plan$wait(-1,c(1,1)),"nonnegative")
fails(plan$wait(0,1),"residual")
fails(plan$wait(0,c(NA,1)),"residual")
fails(plan$wait(0,matrix(1,2,1)),"residual")
fails({plan$step <- identity},"locked")
fails({plan$extra <- TRUE},"locked")
stopifnot(identical(saved,.Random.seed))
cat(comparisons,"prepared/public step comparisons; inversion, captured inputs, independent calls, serialization, validation and RNG contracts passed.\n")
