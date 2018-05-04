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
    dimens <- length(object@mu)
    if (dimens == 1){
      lamb0 <- (object@mu[1]*object@beta[1]/(object@beta[1]- object@alpha[1]) - object@mu[1])/2
      LAMBDA0 <- matrix(rep(lamb0, dimens^2), nrow=dimens, byrow=TRUE)

    } else {
      LAMBDA_st <- solve(diag(dimens) - object@alpha / object@beta) %*% object@mu
      LAMBDA0 <- matrix(rep(LAMBDA_st, dimens), nrow=dimens, byrow=T) * object@alpha / object@beta
    }
    LAMBDA0
  }
)
