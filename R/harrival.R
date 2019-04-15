rarrival <- function(n, mark, type, inter_arrival,
                     N, Nc, lambda, lambda_component,
                     rambda, rambda_component,
                     mu, alpha, beta, dimens, type_col_map, ...){


  # if for each lambda, beta are equal
  test <- 1
  for (i in nrow(beta)){
    test <- test * isTRUE(all.equal(max(beta[i,]), min(beta[i,])))
  }
  if (test == 1){

    #mu can be zero, so warning is turned off for a moment
    oldw <- getOption("warn")
    options(warn = -1)
    candidate_arrival <- stats::rexp(dimens, rate = mu)
    options(warn = oldw)

    current_lambda_component <- matrix(rambda_component[n-1,], nrow = dimens, byrow = TRUE)
    current_lambda <- rowSums(current_lambda_component)

    current_lambda[abs(current_lambda) < 1e-10] <- 0

    matrixD <- 1 + beta[,1] * log(stats::runif(dimens)) / current_lambda
    candidate_arrival <- cbind(candidate_arrival, -1 / beta[,1] * log(pmax(matrixD, 0)))
    inter_arrival <- min(candidate_arrival, na.rm = TRUE)
    minIndex <- which(candidate_arrival == inter_arrival, arr.ind = TRUE)

    type <- minIndex[1]

    if(length(type_col_map) > 0){
      for (j in 1:length(type_col_map)) {
        if (minIndex[1] %in% type_col_map[[j]]){
          type <- j
          break
        }
      }
    }


  } else{

    # thinning algorithm
    current_lambda_component <- matrix(rambda_component[n-1,], nrow = dimens, byrow = TRUE)

    current_lambda <- rowSums(current_lambda_component)
    current_lambda_mu <- current_lambda + mu

    current_lambda_mu[abs(current_lambda_mu) < 1e-10] <- 0


    candidate_arrival <- rep(0, length(current_lambda))

    for(i in 1:length(current_lambda_mu)){

      s <- 0

      if(current_lambda_mu[i] == 0){
        candidate_arrival[i] <- Inf
      } else{
        while(TRUE){

          s <- s + rexp(1, current_lambda_mu[i])

          U <- runif(1, 0, current_lambda_mu[i])

          f <- sum(current_lambda_component[i,] * exp((- beta[i, ] * s)))

          if( U < f + mu[i] ) {
            candidate_arrival[i] <- s
            break
          }

        }
      }
    }


    inter_arrival <- min(candidate_arrival, na.rm = TRUE)
    type <- which(candidate_arrival == inter_arrival)


    # arrival due to components
    # current_lambda_component <- matrix(rambda_component[n-1,], nrow = dimens, byrow = TRUE)
    #
    # matrixD <- 1 + beta * matrix(log(stats::runif(dimens * ncol(beta))), nrow=dimens) / current_lambda_component
    # candidate_arrival <- cbind(candidate_arrival, -1 / beta * log(pmax(matrixD, 0)))
    #
    # print(candidate_arrival)
    #
    # # The minimum is inter arrival time
    # inter_arrival <- min(candidate_arrival, na.rm = TRUE)
    # minIndex <- which(candidate_arrival == inter_arrival, arr.ind = TRUE) #row and col
    #
    # type <- minIndex[1]
    #
    # if(length(type_col_map) > 0){
    #   for (j in 1:length(type_col_map)) {
    #     if (minIndex[1] %in% type_col_map[[j]]){
    #       type <- j
    #       break
    #     }
    #   }
    # }

  }

  c(inter_arrival = inter_arrival, type = type)
}
