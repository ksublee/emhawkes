#' @include hspec.R hmoment.R utilities.R
NULL

#' @export
setGeneric("infer_lambda", function(object, inter_arrival = NULL,
                            type = NULL, mark = NULL,
                            N = NULL, Nc = NULL,
                            lambda0 = NULL, N0 = NULL) standardGeneric("infer_lambda"))
setMethod(
  f="infer_lambda",
  signature(object="hspec"),
  function(object, inter_arrival = NULL,
           type = NULL, mark = NULL,
           N = NULL, Nc = NULL,
           lambda0 = NULL, N0 = NULL){


    #parameter setting
    plist <- setting(object)
    mu <- plist$mu
    alpha <- plist$alpha
    beta <- plist$beta
    impact <- plist$impact
    rmark <- plist$rmark
    dimens <- plist$dimens
    size <-  length(inter_arrival)

    temp <- type_to_N(type, mark, dimens)
    N <- temp[[1]]
    Nc <- temp[[2]]

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

      if (is.function(mu)){
        mu_n <- mu(n = n, mark = mark, type = type, inter_arrival = inter_arrival,
                   N = N, Nc = Nc, lambda = lambda, lambda_component = lambda_component,
                   alpha = alpha, beta = beta)
      } else{
        mu_n <- mu
      }

      # lambda decayed due to time, impact due to mark is not added yet
      decayed <- exp(-beta * inter_arrival[n])
      #decayed_lambda <- current_rambda_component * decayed
      decayed_lambda <- matrix(rambda_component[n-1,], dimens, byrow = T) * decayed

      # update lambda
      lambda_component[n, ] <- t(decayed_lambda)
      lambda[n, ] <- mu_n + rowSums(decayed_lambda)


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

integrate_rambda_component <- function(inter_arrival, rambda_componet, beta, dimens){

  integrated_rambda_component <- matrix(rep(0, nrow(rambda_componet)*ncol(rambda_componet)),
                                            nrow=nrow(rambda_componet))

  current_lambda <- matrix(rambda_componet[1,], nrow = dimens, byrow=TRUE)


  for (n in 2:size) {

    decayed <- exp(-beta * inter_arrival[n])
    integrated_rambda_component[n,] <- t(current_lambda / beta * ( 1 - decayed ))
    new_lambda <- rambda_componet[n,]
    current_lambda <- matrix(new_lambda, nrow=dimens, byrow=TRUE)

  }

  integrated_rambda_component
}


integrate_rambda <- function(inter_arrival, rambda_component, mu, beta, dimens){

  size <- length(inter_arrival)

  if (is.function(mu)){
    # mu is represeted by function
    mu0 <- mu(n = 1, mark = mark, type = type, inter_arrival = inter_arrival,
              N = N, Nc = Nc,
              alpha = alpha, beta = beta)
  } else{
    # mu is a matrix
    mu0 <- mu
  }

  integrated_rambda <- matrix(rep(0, length(rambda_component)/dimens),
                                        nrow=nrow(rambda_component))

  current_lambda <- matrix(rambda_component[1,], nrow = dimens, byrow=TRUE)
  size <- length(inter_arrival)

  for (n in 2:size) {

    decayed <- exp(-beta * inter_arrival[n])

    integrated_rambda[n, ] <- rowSums(current_lambda / beta * ( 1 - decayed )) +
      mu0 * inter_arrival[n]

    new_lambda <- rambda_component[n,]
    current_lambda <- matrix(new_lambda, nrow=dimens, byrow=TRUE)

    # next mu
    if (is.function(mu)){
      # mu is represeted by function
      mu0 <- mu(n = n, mark = mark, type = type, inter_arrival = inter_arrival,
                 N = N, Nc = Nc,
                 alpha = alpha, beta = beta)
    } else{
      # mu is a matrix
      mu0 <- mu
    }

  }

  integrated_rambda

}

#' @export
residual_process <- function(target, type, inter_arrival, rambda_component, mu, beta, dimens){

  integrated_rambda <- integrate_rambda(inter_arrival, rambda_component, mu, beta, dimens)

  #target <- 1
  row_idx <- which(type == target)
  res_process <- rep(0, (length(row_idx)-1))
  for( i in 1:length(res_process)){

    res_process[i] <- sum(integrated_rambda[(row_idx[i] + 1):row_idx[i+1],target])
  }


  res_process

}
