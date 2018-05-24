#' @include hspec.R

setGeneric("get_lambda0", function(object, ...) standardGeneric("get_lambda0"))
#' Get the long-run expectation of intensities - lambda matrix
#'
#' This function is not crucial but useful when you want to set lamba0 by default.
#'
#' @param object hspec
#' @export
setMethod(
  f = "get_lambda0",
  signature(object = "hspec"),
  definition = function(object){


    # parameter setting
    if (is.function(object@mu)){
      mu <- evalf(object@mu)
    } else{
      mu <- object@mu
    }
    if (is.function(object@alpha)){
      alpha <- evalf(object@alpha)
    } else{
      alpha <- object@alpha
    }
    if (is.function(object@beta)){
      beta <- evalf(object@beta)
    } else{
      beta <- object@beta
    }

    dimens <- length(mu)

    if (dimens == 1){
      lamb0 <- (mu[1]*beta[1]/(beta[1]- alpha[1]) - mu[1])/2
      LAMBDA0 <- matrix(rep(lamb0, dimens^2), nrow=dimens, byrow=TRUE)

    } else {
      LAMBDA_st <- solve(diag(dimens) - alpha / beta) %*% mu
      LAMBDA0 <- matrix(rep(LAMBDA_st, dimens), nrow=dimens, byrow=T) * alpha / beta
    }
    LAMBDA0
  }
)
