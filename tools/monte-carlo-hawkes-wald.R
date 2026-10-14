#!/usr/bin/env Rscript
# Run from repository root after R CMD INSTALL .
# Rscript tools/monte-carlo-hawkes-wald.R output-dir reps workers seed
library(emhawkes)
args <- commandArgs(TRUE)
out <- if (length(args)) args[1] else "tools/monte-carlo-results"
reps <- if (length(args) >= 2) as.integer(args[2]) else 200L
workers <- if (length(args) >= 3) as.integer(args[3]) else 1L
seed <- if (length(args) >= 4) as.integer(args[4]) else 20261014L
stopifnot(reps >= 2, workers >= 1, !is.na(seed))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(emhawkes.engine = "auto")
constraints <- list(ineqA = diag(3), ineqB = rep(-1e-6, 3))
# Positivity only: stationarity is a DGP property, not an estimation constraint.
# Generic starts do not use the true excitation or the realized data.
starts <- list(c(.5, .15, 1), c(.5, .7, 1))
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
postprocess <- function(fit) {
  metadata <- attr(fit, "emhawkes")
  ans <- quiet(getFromNamespace("finalize_hawkes_fit", "emhawkes")(
    fit, metadata$specification, metadata$constraints, NULL, list(), "BFGS", "stable"))
  metadata$diagnostics <- ans$diagnostics
  attr(ans$fit, "emhawkes") <- metadata
  ans$fit
}
# The internal helper is used only by this tool to avoid a second optimization.
# Verify the paired public API at one deterministic Hawkes sample.
set.seed(seed)
probe <- hsim(new("hspec", mu=.5, alpha=.5, beta=1), size=101,
              lambda_component0=0)
h <- new("hspec", mu=.5, alpha=.15, beta=1)
a <- quiet(hfit(h, probe$inter_arrival, lambda_component0=0,
                method="BFGS", constraints=constraints, final_hessian="maxLik"))
b <- quiet(hfit(h, probe$inter_arrival, lambda_component0=0,
                method="BFGS", constraints=constraints, final_hessian="stable"))
p <- postprocess(a)
for (field in c("estimate", "maximum", "gradient", "code", "iterations"))
  stopifnot(identical(a[[field]], b[[field]]), identical(a[[field]], p[[field]]))
stopifnot(identical(vcov(b), vcov(p)), identical(fit_diagnostics(b), fit_diagnostics(p)))

settings <- expand.grid(events=c(100L, 500L, 2000L), rho=c(.05, .5, .95))
settings$id <- seq_len(nrow(settings))
settings$mu <- .5; settings$beta <- 1; settings$alpha <- settings$rho
settings$distance <- 1-settings$rho
write.csv(settings, file.path(out, "settings.csv"), row.names=FALSE)
capture.output(sessionInfo(), file=file.path(out,"session-info.txt"))
writeLines(c(paste("seed",seed), paste("replications",reps), paste("workers",workers),
  paste("package",packageVersion("emhawkes")),
  paste("git",paste(system("git rev-parse HEAD",intern=TRUE),collapse=" "))),
  file.path(out,"run-info.txt"))
run_one <- function(r, s) {
  # Each setting/replicate has its own seed, independent of worker scheduling.
  this_seed <- seed + 100000L*s$id + r
  set.seed(this_seed)
  truth <- c(mu=s$mu, alpha=s$alpha, beta=s$beta, rho=s$rho)
  path <- hsim(new("hspec", mu=s$mu, alpha=s$alpha, beta=s$beta),
               size=s$events+1L, lambda_component0=0)
  stopifnot(length(path$inter_arrival)==s$events+1L)
  errors <- character()
  fits <- lapply(starts, function(v) tryCatch(quiet(hfit(
    new("hspec",mu=v[1],alpha=v[2],beta=v[3]), path$inter_arrival,
    lambda_component0=0, method="BFGS", constraints=constraints,
    control=list(reltol=1e-9, iterlim=1000), final_hessian="maxLik")),
    error=function(e) {errors <<- c(errors,conditionMessage(e)); NULL}))
  usable <- vapply(fits,function(f) !is.null(f) && all(is.finite(coef(f))) &&
    is.finite(as.numeric(logLik(f))),logical(1))
  converged <- vapply(fits,function(f) !is.null(f) && maxLik::returnCode(f)==0L,logical(1)) & usable
  pool <- which(if (any(converged)) converged else usable)
  fit <- if (length(pool)) fits[[pool[which.max(vapply(fits[pool],
    function(f) as.numeric(logLik(f)),numeric(1)))]]] else NULL
  opt_ok <- any(converged)
  estimate <- if (is.null(fit)) rep(NA_real_,4) else {
    v <- unname(coef(fit)); c(v, v[2]/v[3])
  }
  names(estimate) <- names(truth)
  stable_error <- ""
  stable <- if (is.null(fit)) NULL else tryCatch(postprocess(fit),
    error=function(e) {stable_error <<- conditionMessage(e); NULL})
  if (!is.null(stable)) for (field in c("estimate","maximum","gradient","code","iterations"))
    stopifnot(identical(fit[[field]],stable[[field]]))
  rows <- lapply(c("maxLik","stable"),function(mode) {
    f <- if (mode=="maxLik") fit else stable
    d <- if (is.null(f)) list() else fit_diagnostics(f)
    V <- if (is.null(f)) NULL else tryCatch(quiet(vcov(f)),error=function(e) NULL)
    # A finite diagonal is insufficient: require positive-definite joint information.
    pd <- !is.null(V) && all(is.finite(V)) &&
      tryCatch(min(eigen((V+t(V))/2,symmetric=TRUE,only.values=TRUE)$values)>0,
               error=function(e) FALSE)
    curvature_ok <- pd && (mode=="maxLik" || identical(d$status,"ok"))
    se <- rep(NA_real_,4)
    if (curvature_ok) {
      se[1:3] <- sqrt(diag(V))
      g <- c(0,1/estimate[3],-estimate[2]/estimate[3]^2)
      se[4] <- sqrt(drop(crossprod(g,V %*% g)))
    }
    covered <- if (opt_ok && curvature_ok) abs(estimate-truth)<=qnorm(.975)*se else rep(NA,4)
    data.frame(setting=s$id, replicate=r, seed=this_seed, method=mode,
      parameter=names(truth), truth=truth, estimate=estimate, se=se, covered=covered,
      optimizer_ok=opt_ok, curvature_ok=curvature_ok,
      optimizer_code=if(is.null(fit)) NA_integer_ else maxLik::returnCode(fit),
      failed_starts=sum(!converged),
      curvature_reason=if(is.null(f)) "fit_unavailable" else if(mode=="maxLik")
        if(pd) "positive_covariance_only" else "invalid_covariance" else
        if(is.null(d$reason)) "ok" else d$reason,
      stationary_estimate=is.finite(estimate[4]) && estimate[4]<1,
      horizon=sum(path$inter_arrival),
      newton_decrement=if(is.null(d$newton_decrement)) NA_real_ else d$newton_decrement,
      error=paste(c(errors,stable_error),collapse=" | "), row.names=NULL)
  })
  do.call(rbind,rows)
}
all_rows <- list()
for (i in seq_len(nrow(settings))) {
  s <- settings[i,]
  t0 <- proc.time()[3]
  result <- if (.Platform$OS.type=="unix" && workers>1) parallel::mclapply(
    seq_len(reps),run_one,s=s,mc.cores=workers,mc.set.seed=FALSE) else
    lapply(seq_len(reps),run_one,s=s)
  stopifnot(all(vapply(result,is.data.frame,logical(1))))
  raw <- do.call(rbind,result)
  write.csv(raw,file.path(out,sprintf("setting-%02d.csv",i)),row.names=FALSE)
  all_rows[[i]] <- raw
  cat(sprintf("setting %d/%d: events=%d rho=%.2f reps=%d %.1fs\n",
      i,nrow(settings),s$events,s$rho,reps,proc.time()[3]-t0))
  flush.console()
}
raw <- do.call(rbind,all_rows)
write.csv(raw,file.path(out,"replicates.csv"),row.names=FALSE)
mean_or_na <- function(x) if(length(x)) mean(x) else NA_real_
summary_rows <- lapply(split(raw,interaction(raw$setting,raw$method,raw$parameter,drop=TRUE)),function(z) {
  ok <- z$optimizer_ok & is.finite(z$estimate)
  valid <- ok & z$curvature_ok & is.finite(z$se)
  common_ids <- with(raw,replicate[setting==z$setting[1] & parameter==z$parameter[1] &
    optimizer_ok & curvature_ok])
  common_ids <- as.integer(names(which(table(common_ids)==2)))
  paired <- valid & z$replicate %in% common_ids
  cover <- z$covered[valid]
  p <- mean_or_na(cover)
  data.frame(setting=z$setting[1],method=z$method[1],parameter=z$parameter[1],
    replications=nrow(z),n_converged=sum(ok),n_valid=sum(valid),n_paired=sum(paired),
    bias=mean_or_na(z$estimate[ok]-z$truth[ok]), empirical_sd=sd(z$estimate[ok]),
    bias_valid=mean_or_na(z$estimate[valid]-z$truth[valid]),
    empirical_sd_valid=sd(z$estimate[valid]),mean_se=mean_or_na(z$se[valid]),
    coverage=p, coverage_mcse=if(length(cover)) sqrt(p*(1-p)/length(cover)) else NA_real_,
    mean_se_paired=mean_or_na(z$se[paired]),coverage_paired=mean_or_na(z$covered[paired]),
    optimizer_failure_rate=mean(!z$optimizer_ok),
    curvature_failure_given_convergence=if(any(ok)) mean(!z$curvature_ok[ok]) else NA_real_,
    interval_available_rate=mean(valid),
    covered_and_available_rate=mean(!is.na(z$covered) & z$covered),
    nonstationary_given_convergence=if(any(ok)) mean(!z$stationary_estimate[ok]) else NA_real_)
})
summary <- merge(do.call(rbind,summary_rows),settings,by.x="setting",by.y="id")
write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE)
reasons <- as.data.frame(with(raw[raw$parameter=="mu",],
  table(setting,method,optimizer_ok,curvature_reason)))
write.csv(reasons,file.path(out,"failure-reasons.csv"),row.names=FALSE)
