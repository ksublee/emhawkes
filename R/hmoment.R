#' @include hspec.R

setGeneric("get_lambda0", function(object, ...) standardGeneric("get_lambda0"))
# Get the long-run expectation of intensities - lambda matrix
#
# This function is not crucial but useful when you want to set lamba0 by default.
#
# @param object hspec
setMethod(
  f = "get_lambda0",
  signature(object = "hspec"),
  definition = function(object, mark = mark, type = type, inter_arrival = inter_arrival,
                        N = N, Nc = Nc,
                        mu = mu, alpha = alpha, beta = beta,
                        lambda = NULL, lambda_component = NULL, lambda_component_n = NULL,
                        ...){

    dimens <- object@dimens

    # parameter setting
    if (is.function(object@mu)){
      if(length(formals(object@mu)) == 1){
        mu0 <- evalf(object@mu)
      } else {
        mu0 <- mu(n = 1, mark = mark, type = type, inter_arrival = inter_arrival,
                  N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
                  lambda_component_n = lambda_component_n,
                  alpha = alpha, beta = beta)
      }
    } else{
      mu0 <- object@mu
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


    if (dimens == 1){
      lamb0 <- (mu[1]*beta[1]/(beta[1]- alpha[1]) - mu[1])/2
      LAMBDA0 <- matrix(lamb0, nrow=1)

    } else {
      # need to handle when the matrix is singular
      LAMBDA0 <- matrix(rep(0, dimens^2), nrow=dimens)
      LAMBDA0 <- tryCatch({
          LAMBDA_st <- solve(diag(dimens) - alpha / beta) %*% mu0
          matrix(rep(LAMBDA_st, dimens), nrow=dimens, byrow=T) * alpha / beta
        },
        error = function(e){
          warning("Due to the singualr martrix in caculcation, set initial value of lambda to mu.")
          matrix(rep(0, dimens^2), nrow=dimens)
        }
      )
    }
    LAMBDA0
  }
)
