#' @include hspec.R hmoment.R
#'
#' Compute the loglikelihood function
#'
#' This is a generic function.
#' The loglikelihood of the ground process of the Hawkes model.
#' (The estimation for jump distribution is not provided.)
#'
#' @param object \code{\link{hspec-class}}. The parameter values in the object are used to compute the log-likelihood.
#' @param inter_arrival Inter-arrival times of events. Includes inter-arrival for events that occur in all dimensions. Start with zero.
#' @param event_idx a vector of dimensions. Distinguished by numbers, 1, 2, 3, and so on. Start with zero.
#' @param mark a vector of mark (jump) sizes. Start with zero.
#' @param lambda0 The starting values of lambda. Must have the same dimensional matrix (n by n) with \code{hspec}.
#'
#' @seealso \code{\link{hspec-class}}, \code{\link{hfit,hspec-method}}
#' @export
setMethod(
  f="logLik",
  signature(object="hspec"),
  function(object, inter_arrival, event_idx=NULL, mark=NULL, lambda0=NULL){

    # When the mark sizes are not provided, the jumps are all unit jumps
    if(is.null(mark)) {
      mark <- rep(1, length(inter_arrival))
    }

    # dimension of Hawkes process
    dimens <- length(object@mu)

    # if dimens == 1 and event_idx is not provided, then all event_idx is 1.
    if(dimens==1 & is.null(event_idx)) {
      event_idx <- rep(1, length(inter_arrival))
    } else if (dimens != 1 & is.null(event_idx)) {
      stop("The argument event_idx should be provided.")
    }

    # parameter setting

    mu <- object@mu
    alpha <- object@alpha
    beta <- object@beta
    impact <- object@impact

    # default lambda0
    if(is.null(lambda0)) {
      warning("The initial values for intensity processes are not provided. Internally determined initial values are used.\n")
      lambda0 <- get_lambda0(object)
    }

    # n is length(inter_arrival) - 1
    n <- length(inter_arrival)

    #if (dimens==1) rowSums_lambda0 <- lambda0
    #else rowSums_lambda0 <- rowSums(lambda0)
    rowSums_lambda0 <- rowSums(matrix(lambda0, nrow=dimens))

    sum_log_lambda <- 0
    sum_integrated_lambda_component <- 0

    current_lambda <- lambda0

    for (i in 2:n) {


      # update lambda
      # if (dimens == 1) {
      #   Impact <- alpha * (1 + (mark[i+1] - 1) * ETA )
      # } else {
      #   Impact <- matrix(rep(0, dimens^2), nrow = dimens)
      #   Impact[ , event_idx[i+1]] <- alpha[ , event_idx[i+1]] * (1 + (mark[i+1] - 1) * ETA[ , event_idx[i+1]])
      # }
      #


      decayed <- exp(-beta * inter_arrival[i])
      decayed_lambda <- current_lambda * decayed

      # impact by alpha
      impact_alpha <- matrix(rep(0, dimens^2), nrow = dimens)
      impact_alpha[ , event_idx[i]] <- alpha[ , event_idx[i]]

      # new_lambda = [[lambda11, lambda12, ...], [lambda21, lambda22, ...], ...]
      if(!is.null(impact)){
        # impact by mark
        impact_mark <- matrix(rep(0, dimens^2), nrow = dimens)

        impact_res <- impact(i = i, N = N, Ng = Ng, mark = mark,
                             lambda = lambda, lambda_component = lambda_component,
                             realized_mark_lambda = realized_mark_lambda,
                             event_idx = event_idx)
        #impact_res <- impact()

        impact_mark[ , event_idx[i]] <- impact_res[ , event_idx[i]]

        new_lambda <- decayed_lambda + impact_alpha + impact_mark

      } else {

        new_lambda <- decayed_lambda + impact_alpha

      }

      # sum of integrated_lambda_component
      sum_integrated_lambda_component <- sum_integrated_lambda_component + sum(current_lambda / beta * ( 1 - decayed ))

      # sum of log lambda when jump occurs
      if (dimens == 1) lambda_lc <- mu + decayed_lambda
      else lambda_lc <- mu + rowSums(decayed_lambda)
      sum_log_lambda <- sum_log_lambda + log(lambda_lc[event_idx[i]])

      # current_lambda <- matrix(lambda_component[i, ], nrow = dimens, byrow = TRUE)
      current_lambda <- new_lambda  # lambda determined in the previous loop

    }

    # log likelihood for ground process
    sum_log_lambda - sum(mu*sum(inter_arrival)) - sum_integrated_lambda_component

  }
)
