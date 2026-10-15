#' @include hspec.R hmoment.R hinterface.R
NULL

#' Compute the Log-Likelihood Function
#'
#' Evaluate a specified model's likelihood, with optional intensity and residual paths.
#'
#' @param object An \code{\link{hspec-class}} object containing parameter values for computing the log-likelihood.
#' @param inter_arrival A vector of inter-arrival times for events across all dimensions, starting with zero.
#' @param type A vector indicating the dimensions, represented by numbers (1, 2, 3, etc.), starting with zero.
#' @param mark A vector of mark (jump) sizes, starting with zero.
#' @param N Optional cumulative event-count matrix, with one row per input row
#' and one column per event type. Each row contains counts after that event.
#' @param Nc Optional cumulative mark-total matrix with the same shape as `N`.
#' @param lambda_component0 Initial values for the lambda component \eqn{\lambda_{ij}}.
#' A finite scalar is expanded to every kernel component. A matrix must match
#' the shape of `alpha` and `beta`. Supply this argument when `type_col_map` is set.
#' @param N0 Initial event counts: a finite scalar or one value per event type.
#' If omitted, use the first row of supplied `N`, or zero when `N` is absent.
#' @param Nc0 Initial mark totals: a finite scalar or one value per event type.
#' If omitted, use the first row of supplied `Nc`, or zero when `Nc` is absent.
#' @param infer Logical; if `TRUE`, return intensity paths and residual processes
#' together with the likelihood.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return With `infer = FALSE`, a length-one numeric vector named `loglikelihood`.
#' With `infer = TRUE`, a list containing:
#' * `loglikelihood`: the same likelihood value.
#' * `lambda` and `rambda`: intensity matrices immediately before and after events.
#' * `lambda_component` and `rambda_component`: the corresponding kernel components,
#'   with columns ordered by intensity row, then kernel column.
#' * `res_process1`, ..., `res_processD`: residuals between successive events of
#'   each type. Types with fewer than two events have empty residual vectors.
#'
#' Intensity matrices include the initial row. For residual-distribution models,
#' `res_process` values use the specified distribution and `eres_process1`, ...,
#' `eres_processD` give their unit-exponential transforms.
#'
#' @details
#' The first row represents the initial state and is not an event.
#' Omit `type` for a univariate model and `mark` for unit marks. The observation
#' vectors have the same length. The likelihood ends at the last supplied event.
#' Counting-process rows contain values after the corresponding event.
#' Supplied `N` and `Nc` matrices must agree with `type`, `mark`, and any
#' supplied initial counts. Event order is preserved, including zero intervals.
#' Component functions receive the arguments they declare; functions with `...`
#' receive the available event context. Intensity histories are `NULL` before
#' initialization. Baseline functions use `n - 1` to access the preceding event.
#' If initial intensity components are omitted, they are determined from the
#' supplied model. For a fitted model, use [fitted_hspec()] with the arguments
#' returned by [fit_initial_state()] to preserve the fit's initialization setting.
#'
#' With the default automatic selection, standalone `logLik()` uses the Rcpp
#' engine for numeric Hawkes specifications. Function-valued parameters,
#' impact functions, mark densities, residual distributions, and classed numeric
#' inputs use the R engine. `rmark` and `rresidual` are simulation generators and
#' do not by themselves prevent numeric likelihood acceleration.
#' Set `options(emhawkes.engine = "r")` for the R calculation, `"auto"` for
#' automatic selection, or `"cpp"` to require a supported numeric specification.
#' Setting `"cpp"` is not needed for automatic acceleration. The engine supports
#' unequal decay rates within a row and multi-kernel mappings. During fitting,
#' [hfit()] also supports C++ adapters for the callback models documented there.
#' Engine selection does not change simulation or enable analytic gradients.
#'
#' @seealso [fitted_hspec()], [fit_initial_state()],
#' \code{\link{hspec-class}}, \code{\link{hfit,hspec-method}}
#'
#' @docType methods
#' @rdname logLik
#' @aliases logLik
#'
#' @export
setMethod(
  f="logLik",
  signature(object="hspec"),
  function(object, inter_arrival, type = NULL, mark = NULL, N = NULL, Nc = NULL,
           N0 = NULL, Nc0 = NULL, lambda_component0 = NULL,
           infer = FALSE, ...){

    additional_argument <- list(...)
    if ("lambda0" %in% names(additional_argument)) {

      warning("lambda0 is deprecated; instead use lambda_component0.")

      lambda_component0 <- additional_argument[["lambda0"]]

    }

    callbacks <- model_callback_plans(object)
    observation <- normalize_observation(object, inter_arrival, type, mark, N, Nc, N0, Nc0)
    loglik_hspec(object, observation, lambda_component0, infer, callbacks,
                 show_initial_message = !("showWarning" %in% names(additional_argument)))
  }
)

# Internal entry point for observations already checked by hfit.
loglik_hspec <- function(object, observation, lambda_component0 = NULL, infer = FALSE,
                         callbacks = model_callback_plans(object), show_initial_message = TRUE,
                         engine = getOption("emhawkes.engine", "auto"), execution = NULL) {
    engine <- match.arg(engine, c("auto", "r", "cpp"))
    if (engine != "r") {
      reason <- native_hawkes_reason(object, observation)
      if (is.null(reason)) {
        record_fit_engine(execution, "cpp")
        return(loglik_hspec_cpp(object, observation, lambda_component0, infer,
                                show_initial_message))
      }
      if (engine == "cpp") stop("C++ Hawkes engine does not support this model: ", reason, call. = FALSE)
    }
    record_fit_engine(execution, "r", if (engine == "auto") reason else NULL)
    loglik_hspec_r(object, observation, lambda_component0, infer, callbacks,
                    show_initial_message)
}

# Keep the R implementation as the general callback engine and reference.
loglik_hspec_r <- function(object, observation, lambda_component0 = NULL, infer = FALSE,
                           callbacks = model_callback_plans(object), show_initial_message = TRUE) {

    ## Parameter setting and initialization
    # after parameter setting, functions mu, alpha, beta become matrices, if possible
    # Note that mu, alpha, beta, eta can be given functions even they are essentially matrices.
    plist <- setting(object)
    # mu, alpha, beta, eta
    list2env(plist, envir = environment())

    inter_arrival <- observation$inter_arrival
    type <- observation$type
    mark <- observation$mark
    N <- observation$N
    Nc <- observation$Nc
    size <- length(inter_arrival)
    lambda <- lambda_component <- lambda_component_n <- NULL
    mu_call <- prepare_component_call(callbacks$mu)
    impact_call <- prepare_component_call(callbacks$impact)
    dmark_call <- prepare_component_call(callbacks$dmark)

    # check the argument lists in mu and impact
    mu_args <- impct_args <- c()
    if(is.function(mu)) mu_args <- methods::formalArgs(mu)
    if(!is.null(object@impact)) impct_args <- methods::formalArgs(object@impact)
    mark_args <- if (is.function(object@dmark)) methods::formalArgs(object@dmark) else character(0)
    history_args <- c(impct_args, mu_args, mark_args)
    impact_has_dots <- "..." %in% impct_args


    if (is.function(mu)){

      mu0 <- mu_call(mu, 1, mark, type, inter_arrival, N, Nc, lambda,
                   lambda_component, lambda_component_n, mu, alpha, beta)


    } else{
      # mu is a matrix
      mu0 <- mu
    }

    # ---- Lambda component0 check & setup ----
    if (is.null(lambda_component0)) {
      if (length(object@type_col_map) > 0) {
        stop("In this model, please provide lambda_component0.")
      }

      if (show_initial_message) {
        message("The initial values for intensity processes are not provided. Internally determined initial values are used.\n")
      }

      lambda_component0 <- get_lambda0(
        object,
        mark = mark,
        type = type,
        inter_arrival = inter_arrival,
        N = N, Nc = Nc,
        mu = mu,
        alpha = alpha,
        beta = beta
      )
    }
    lambda_component0 <- normalize_lambda_initial(lambda_component0, object@dimens, ncol(beta))
    lambda_component_n <- lambda_component0

    # ---- Determine necessity for lambda/lambda_component ----
    is_lambda_necessary <- infer || "lambda" %in% history_args
    is_lambda_component_necessary <- infer || "lambda_component" %in% history_args

    # ---- lambda_component ----
    if (is_lambda_component_necessary) {
      lambda_component <- initialize_lambda_matrix(
        prefix = "lambda",
        size = size,
        dimens = object@dimens,
        beta_cols = ncol(beta),
        lambda0 = lambda_component0
      )
    }

    # ---- lambda ----
    if (is_lambda_necessary) {
      lambda <- initialize_lambda_total(
        prefix = "lambda",
        size = size,
        dimens = object@dimens,
        mu0 = mu0,
        lambda0 = lambda_component0
      )
    }

    # ---- rambda (inference mode) ----
    if (infer) {
      rambda_component <- initialize_lambda_matrix(
        prefix = "rambda",
        size = size,
        dimens = object@dimens,
        beta_cols = ncol(beta),
        lambda0 = lambda_component0
      )

      rambda <- initialize_lambda_total(
        prefix = "rambda",
        size = size,
        dimens = object@dimens,
        mu0 = mu0,
        lambda0 = lambda_component0
      )

      # phi matrix (phi(tau) = epsilon), so this is related to residual process
      phis <- matrix(nrow = length(inter_arrival), ncol = object@dimens)
      colnames(phis) <- paste0("phi", 1:object@dimens)

      # historical right-continuous version of mus
      rmus <- matrix(nrow = length(inter_arrival), ncol = object@dimens)
      rmus[1, ] <- mu0
    }


    # ---- Initialize accumulators ----
    accumulators <- c("sum_log_lambda", "sum_integrated_lambda_component",
                      "sum_mu_inter_arrival", "sum_log_f_epsilon_phi",
                      "sum_log_dmark")

    for (name in accumulators) assign(name, 0)


    # ---- Evaluate piecewise constant mu (only supported form) ----
    rmu_n <- if (is.function(mu) && size > 1L) {
      mu_call(mu, 2, mark, type, inter_arrival, N, Nc, lambda,
                   lambda_component, lambda_component_n, mu, alpha, beta)
    } else {
      mu0
    }


    # ---- Preallocate zero matrices for alpha, eta, and impact ----
    zero_mat_for_alpha <- init_zero_matrix(object@dimens, ncol(beta))
    zero_mat_for_eta <- zero_mat_for_alpha  # same size
    zero_mat_for_impact <- zero_mat_for_alpha  # same size

    rambda_component_n <- lambda_component0

    # Decay and integration weights depend only on beta and observed intervals.
    # Bound the temporary tables; custom numeric classes keep their own operators.
    cache_kernels <- !is.object(beta) && !is.object(inter_arrival) &&
      length(beta) <= 8192L && size >= 32L
    cache_end <- 1L
    if (cache_kernels) {
      kernel_attributes <- attributes(beta)
      scalar_kernel <- length(beta) == 1L && identical(names(kernel_attributes), "dim")
      rates <- as.vector(beta)
      block_size <- max(1L, 65536L %/% length(rates))
    }

    ## main loop
    for (n in seq_len(size - 1L) + 1L) {

      mu_n <- rmu_n; rambda_component_n_minus_1 <- rambda_component_n

      type_n <- type[n]

      if (cache_kernels) {
        if (n > cache_end) {
          cache_start <- n - 1L
          cache_end <- n + min(block_size - 1, size - n)
          block_dt <- inter_arrival[seq.int(n, cache_end)]
          scaled_time <- outer(rates, block_dt)
          cached_decay <- exp(-scaled_time)
          cached_weight <- -expm1(-scaled_time) / rates
          zero <- which(scaled_time == 0)
          cached_weight[zero] <- block_dt[(zero - 1L) %/% length(rates) + 1L]
        }
        index <- n - cache_start
        decayed <- cached_decay[, index]
        weight <- cached_weight[, index]
        if (!scalar_kernel) {
          attributes(decayed) <- kernel_attributes
          attributes(weight) <- kernel_attributes
        }
        lambda_component_n <- rambda_component_n_minus_1 * decayed
        integrated_component <- rambda_component_n_minus_1 * weight
        integrated_component[rambda_component_n_minus_1 == 0] <- 0
      } else {
        decayed <- exp(-beta * inter_arrival[n])
        lambda_component_n <- rambda_component_n_minus_1 * decayed
        integrated_component <- integrate_exp_kernel(rambda_component_n_minus_1,
                                                      beta, inter_arrival[n])
      }

      # integrate mu over inter_arrival w.r.t. time
      if ((is.matrix(mu_n) && nrow(mu_n) == object@dimens && ncol(mu_n) == 1) ||
          (is.vector(mu_n) && length(mu_n) == 1)) {

        # col num : dimens
        integrated_mu <- mu_n * inter_arrival[n]

      } else {

        # When mu is observed as piecewise constant between n-1 and n
        # Rows from 1 to object@dimen are mu values for each dimension.
        # The last row represents the inter-arrival where mu changes between n-1 and n.
        sum_mu_n <- colSums(mu_n[1:object@dimens, ])

        # col num : dimens
        integrated_mu <- sum_mu_n * mu_n[object@dimens + 1, ]

      }

      if (is.null(object@dresidual) && !infer){

        # basic model and do not need to infer
        ## 1. sum of integrated_lambda_component and sum of mu * inter_arrival
        ## In the log-likelihood calculation, a negative sign is applied, so it's subtracted here.
        sum_log_f_epsilon_phi <- sum_log_f_epsilon_phi - sum(integrated_component) - sum(integrated_mu)


      } else if (is.null(object@dresidual)) {

        # basic model and need to infer
        for (k in 1:object@dimens){

          # Integrate each kernel with its own decay rate.
          phi_k_wo_mu <- sum(integrated_component[k, ])

          # partial residual for k-th dimension
          phi_value_k <- phi_k_wo_mu + integrated_mu[k]
          sum_log_f_epsilon_phi <- sum_log_f_epsilon_phi - phi_value_k

          phis[n, k] <- phi_value_k

        }

      } else {


        # flexible model

        for (k in 1:object@dimens){

          # Integrate each kernel with its own decay rate.
          phi_k_wo_mu <- sum(integrated_component[k, ])

          # partial residual for k-th dimension
          phi_value_k <- phi_k_wo_mu + integrated_mu[k]

          log_value <- suppressWarnings(if (k == type_n) {
            residual_density(object@dresidual, phi_value_k, log = TRUE)
          } else {
            residual_cdf(object@presidual, phi_value_k, lower.tail = FALSE, log.p = TRUE)
          })
          sum_log_f_epsilon_phi <- sum_log_f_epsilon_phi + log_value

          if(infer) phis[n, k] <- phi_value_k

        }


      }



      ## 2. sum of log lambda when jump occurs
      # Calculate sum of log lambda when a jump occurs
      lambda_at_event <- mu_n[type_n] + sum(lambda_component_n[type_n, ])

      # Nonnegative ordinary intensities cannot produce a log-domain warning.
      # Retain warning suppression for invalid proposals and custom Math methods.
      log_lambda <- if (!is.object(lambda_at_event) && is.double(lambda_at_event) &&
                        isTRUE(lambda_at_event >= 0)) {
        log(lambda_at_event)
      } else suppressWarnings(log(lambda_at_event))
      sum_log_lambda <- sum_log_lambda + log_lambda


      if(is_lambda_component_necessary) lambda_component[n, ] <- t(lambda_component_n)
      if(is_lambda_necessary) lambda[n, ] <- as.vector(mu_n) + rowSums(lambda_component_n)

      ## 3. sum of log mark p.d.f. when jump occurs
      if(!is.null(object@dmark)){
        sum_log_dmark <- sum_log_dmark +
          log(dmark_call(object@dmark, n, mark, type, inter_arrival, N, Nc, lambda,
                   lambda_component, lambda_component_n, mu, alpha, beta))

      }



      ## impact

      # Determine types based on the type_col_map
      types <- if (is.null(object@type_col_map)) {
        type_n
      } else if (length(object@type_col_map) > 0) {
        object@type_col_map[[type_n]]
      } else {
        stop("Check the type_col_map argument.")
      }

      # 1. impact by alpha

      impact_alpha <- zero_mat_for_alpha   # matrix
      impact_alpha[ , types] <- alpha[ , types]

      rambda_component_n <- lambda_component_n + impact_alpha

      # 2. additional impact by eta

      if(!is.null(eta)) {

        impact_eta <- zero_mat_for_eta   # matrix
        impact_eta[ , types] <- eta[ , types] * (mark[n] - 1)

        rambda_component_n <- rambda_component_n + impact_eta
      }


      # 3. impact by impact function
      # new_lambda = [[lambda11, lambda12, ...], [lambda21, lambda22, ...], ...]
      if(!is.null(object@impact)){

        impact_mark <- zero_mat_for_impact
        impact_res <- if (impact_has_dots) {
          forceAndCall(12L, object@impact, n = n, mark = mark, type = type,
                       inter_arrival = inter_arrival, N = N, Nc = Nc, lambda = lambda,
                       lambda_component = lambda_component, lambda_component_n = lambda_component_n,
                       mu = mu, alpha = alpha, beta = beta)
        } else impact_call(object@impact, n, mark, type, inter_arrival, N, Nc, lambda,
                           lambda_component, lambda_component_n, mu, alpha, beta)

        impact_mark[ , types] <- impact_res[ , types]

        # for next step
        rambda_component_n <- rambda_component_n + impact_mark

      }

      # new mu, i.e., right continuous version of mu
      if (is.function(mu)){
        # mu is represented by a function
        rmu_n <- mu_call(mu, n + 1, mark, type, inter_arrival, N, Nc, lambda,
                   lambda_component, lambda_component_n, mu, alpha, beta)
      } else{
        # mu is a matrix
        rmu_n <- mu
      }

      if(infer){
        rambda_component[n, ] <- t(rambda_component_n)
        rambda[n, ] <- as.vector(rmu_n) + rowSums(rambda_component_n)
        rmus[n, ] <- rmu_n
      }



    }


    if(!infer) {
      # only log-likelihood
      c(loglikelihood = sum_log_lambda + sum_log_f_epsilon_phi + sum_log_dmark)

    } else {

      # residual inference

      return_object <- list(loglikelihood = sum_log_lambda + sum_log_f_epsilon_phi + sum_log_dmark,
                            lambda = lambda, lambda_component = lambda_component,
                            rambda = rambda, rambda_component = rambda_component)

      if (is.null(object@dresidual)){

        # basic model

        for (k in 1:object@dimens){

          # the ordered sequence of indices associated with type k
          n_k <- which(type == k & seq_along(type) > 1L)
          res_process <- numeric(max(0L, length(n_k) - 1L))

          for(j in seq_along(res_process)){
            # n_k[j]+1 : index of previous jump + 1
            # n_k[j+1] : index of current jump
            res_process[j] <- sum(phis[(n_k[j]+1):n_k[j+1], k])
          }

          return_object[[paste0("res_process", k)]] <- res_process

        }

      } else {

        # flexible model
        for (k in 1:object@dimens){

          # the ordered sequence of indices associated with type k
          n_k <- which(type == k & seq_along(type) > 1L)
          res_process <- numeric(max(0L, length(n_k) - 1L))
          eres_process <- numeric(max(0L, length(n_k) - 1L))

          for(j in seq_along(res_process)){


            # n_k[j]+1 : index of previous jump + 1
            # n_k[j+1] : index of current jump
            # F(phi(tau)) = F(epsilon)
            log_survival <- sum(residual_cdf(object@presidual,
              phis[(n_k[j]+1):n_k[j+1], k], lower.tail = FALSE, log.p = TRUE))
            res_process[j] <- residual_quantile(object@qresidual, log_survival,
                                               lower.tail = FALSE, log.p = TRUE)
            eres_process[j] <- -log_survival

          }

          return_object[[paste0("res_process", k)]] <- res_process
          return_object[[paste0("eres_process", k)]] <- eres_process

        }
      }

      return_object

    }

}

# Helpers: Initialization for lambda matrices
initialize_lambda_matrix <- function(prefix, size, dimens, beta_cols, lambda0) {
  mat <- matrix(0, nrow = size, ncol = dimens * beta_cols)
  mat[1, ] <- t(lambda0)
  colnames(mat) <- paste0(prefix, rep(1:dimens, each = beta_cols), 1:beta_cols)
  return(mat)
}

initialize_lambda_total <- function(prefix, size, dimens, mu0, lambda0) {
  mat <- matrix(0, nrow = size, ncol = dimens)
  mat[1, ] <- as.vector(mu0) + rowSums(matrix(lambda0, nrow = dimens))
  colnames(mat) <- paste0(prefix, 1:dimens)
  return(mat)
}

# ---- Preallocate zero matrices for alpha, eta, and impact ----
init_zero_matrix <- function(dimens, beta_cols) {
  matrix(0, nrow = dimens, ncol = beta_cols)
}

