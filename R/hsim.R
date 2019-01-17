#' @include hspec.R hmoment.R utilities.R
NULL

#' Simulate a multivariate Hawkes process.
#' Generic function hsim.
#'
#'
#' @param object \code{\link{hspec-class}}. This object includes the parameter values.
#' @param lambda0 the starting values of lambda component. numeric or matrix.
#' @param N0 the starting values of N
#' @param size the number of observations.
#'
#' @return hreal S3-object, summary of the realization of the Haweks model
#'
#' @rdname hsim
#'
#' @examples
#' mu <- c(0.1, 0.1)
#' alpha <- matrix(c(0.2, 0.1, 0.1, 0.2), nrow=2, byrow=TRUE)
#' beta <- matrix(c(0.9, 0.9, 0.9, 0.9), nrow=2, byrow=TRUE)
#' h <- new("hspec", mu=mu, alpha=alpha, beta=beta)
#' res <- hsim(h, size=100)
#'
#' @export
setGeneric("hsim", function(object, size = 100, lambda0 = NULL, N0 = NULL) standardGeneric("hsim"))
#'
#' The method simulate multivariate Hawkes processes.
#' The object \code{\link{hspec-class}} contains the parameter values such as \code{mu}, \code{alpha}, \code{beta}.
#' The mark (jump) structure may or may not be included.
#' It returns an object of class \code{hreal} which contains \code{inter_arrival}, \code{arrival},
#' \code{type}, \code{mark}, \code{N}, \code{Nc}, \code{lambda}, \code{lambda_component}, \code{rambda}, \code{rambda_component}.
#'
#' @rdname hsim
setMethod(
  f="hsim",
  signature(object = "hspec"),
  definition = function(object, size = 100, lambda0 = NULL, N0 = NULL){

    # parameter setting
    plist <- setting(object)
    mu <- plist$mu
    alpha <- plist$alpha
    beta <- plist$beta
    impact <- plist$impact
    rmark <- plist$rmark
    dimens <- plist$dimens

    N <- matrix(numeric(length = dimens * size), ncol = dimens)
    Nc  <- matrix(numeric(length = dimens * size), ncol = dimens)
    colnames(Nc)  <- paste0("Nc", 1:dimens)
    colnames(N) <- paste0("N", 1:dimens)

    type <- numeric(length = size)
    inter_arrival <- numeric(length = size)
    mark <- numeric(length = size)

    if (!is.null(N0)){
      N[1,] <- N0
      Nc[1,] <- N0
    }

    if (is.function(mu)){
      mu0 <- mu(n = 1, mark = mark, type = type, inter_arrival = inter_arrival,
                 N = N, Nc = Nc,
                 alpha = alpha, beta = beta)
    } else {
      mu0 <- mu
    }

    if(!is.null(lambda0)){

      # If the dimensions of model and lambda0 do not match, lambda0 will be adjusted
      if (dimens^2 > length(lambda0)){
        warning("The size of lambda0 does not match to the dimension of the model and is adjusted. \n
                lambda0 is now :")
        lambda0 <- rep(lambda0, dimens^2)
      }
      if (dimens^2 < length(lambda0)){
        warning("The size of lambda0 does not match to the dimension of the model and is adjusted.\n
                lambda0 is now :")
        lambda0 <- lambda0[1:dimens^2]
      }

      lambda0 <- as.matrix(lambda0, nrow = dimens)

    } else {
      # default lambda0
      warning("The initial values for intensity processes are not provided. Internally determined initial values are used.\n")
      lambda0 <- get_lambda0(object, mark = mark, type = type, inter_arrival = inter_arrival,
                             N = N, Nc = Nc,
                             alpha = alpha, beta = beta)
    }

    # Preallocation for lambdas and Ns and set initial values for lambdas
    lambda_component <- matrix(sapply(lambda0, c, numeric(length = size - 1)), ncol = dimens^2)
    rowSums_lambda0 <- rowSums(matrix(lambda0, nrow=dimens))

    lambda <- matrix(sapply(mu0 + rowSums_lambda0, c, numeric(length = size - 1)), ncol = dimens)

    rambda <- lambda
    rambda_component <- lambda_component

    # Set column names
    colnames(lambda) <- paste0("lambda", 1:dimens)
    indxM <- matrix(rep(1:dimens, dimens), byrow = TRUE, nrow = dimens)
    colnames(lambda_component) <- paste0("lambda", indxM, t(indxM))

    colnames(rambda) <- paste0("rambda", 1:dimens)
    indxM <- matrix(rep(1:dimens, dimens), byrow = TRUE, nrow = dimens)
    colnames(rambda_component) <- paste0("rambda", indxM, t(indxM))


    #current_rambda_component <- lambda0
    for (n in 2:size) {

      # only piecewise constant mu is available
      if (is.function(mu)){
        mu_n <- mu(n = n - 1, mark = mark, type = type, inter_arrival = inter_arrival,
                   N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
                   alpha = alpha, beta = beta)
      } else{
        mu_n <- mu
      }

      ### Determine the next arrival #############################################################
      res <- rarrival(n = n, mark = mark, type = type, inter_arrival = inter_arrival,
               N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
               rambda = rambda, rambda_component = rambda_component,
               mu = mu_n, alpha = alpha, beta = beta, dimens = dimens)
      inter_arrival[n] <- res["inter_arrival"]
      type[n] <- res["type"]
      ############################################################################################

      N[n, ] <- N[n-1, ]
      N[n, type[n]] <- N[n-1, type[n]] + 1

      # lambda decayed due to time, impact due to mark is not added yet
      decayed <- exp(-beta * inter_arrival[n])
      #decayed_lambda <- current_rambda_component * decayed
      decayed_lambda <- matrix(rambda_component[n-1,], dimens, byrow = T) * decayed

      # update lambda
      lambda_component[n, ] <- t(decayed_lambda)
      lambda[n, ] <- mu_n + rowSums(decayed_lambda)



      # generate a mark for Hawkes
      # This quantity is added to the counting process.
      if( !is.null(rmark) ){
        # mark may depends on other variables
        # mark[n] is a scalar
        mark[n] <- rmark(n = n, Nc = Nc, N = N,
                               lambda = lambda, lambda_component = lambda_component,
                               type = type)
      } else {
        mark[n] <- 1
      }

      Nc[n, ] <- Nc[n-1, ]
      Nc[n, type[n]] <- Nc[n-1, type[n]] + mark[n]


      # impact by alpha
      impact_alpha <- matrix(rep(0, dimens^2), nrow = dimens)
      impact_alpha[ , type[n]] <- alpha[ , type[n]]

      # new_lambda = [[lambda11, lambda12, ...], [lambda21, lambda22, ...], ...]
      if(!is.null(impact)){
        # impact by mark
        impact_mark <- matrix(rep(0, dimens^2), nrow = dimens)

        impact_res <- impact(n = n, mark = mark, type = type, inter_arrival = inter_arrival,
                             N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
                             mu = mu, alpha = alpha, beta = beta)

        impact_mark[ , type[n]] <- impact_res[ , type[n]]

        new_lambda <- decayed_lambda + impact_alpha + impact_mark

      } else {

        new_lambda <- decayed_lambda + impact_alpha

      }

      # update rambda
      # rambda_component = {"rambda11", "rambda12", ..., "rambda21", "rambda22", ...}
      rambda_component[n, ] <- t(new_lambda)

      if (is.function(mu)){
        mu_n <- mu(n = n, mark = mark, type = type, inter_arrival = inter_arrival,
                   N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
                   alpha = alpha, beta = beta)
      } else{
        mu_n <- mu
      }
      rambda[n, ] <- mu_n + rowSums(new_lambda)

      #current_rambda_component <- new_lambda
    }

    realization <- list(object, inter_arrival, cumsum(inter_arrival), type, mark,
                        N, Nc, lambda, lambda_component, rambda, rambda_component)
    names(realization) <- c("hspec", "inter_arrival", "arrival", "type", "mark",
                            "N", "Nc", "lambda", "lambda_component", "rambda", "rambda_component")
    class(realization) <- c("hreal")

    realization
  }
)
