rarrival <- function(n, mark, type, inter_arrival,
                     N, Nc, lambda, lambda_component,
                     rambda, rambda_component,
                     mu, alpha, beta, dimens, ...){
  #mu can be zero, so warning is turned off for a moment
  oldw <- getOption("warn")
  options(warn = -1)
  candidate_arrival <- stats::rexp(dimens, rate = mu)
  options(warn = oldw)
  #current_lambda_component <- matrix(as.numeric(lambda_component[n-1, ]), nrow = dimens, byrow = TRUE)

  # arrival due to components

  current_lambda_component <- matrix(rambda_component[n-1,], nrow = dimens, byrow = TRUE)

  matrixD <- 1 + beta * log(stats::runif(dimens^2)) / current_lambda_component
  candidate_arrival <- cbind(candidate_arrival, -1 / beta * log(pmax(matrixD, 0)))

  # The minimum is inter arrival time
  inter_arrival <- min(candidate_arrival, na.rm = TRUE)
  minIndex <- which(candidate_arrival == inter_arrival, arr.ind = TRUE) #row and col

  type <- minIndex[1]  # row

  c(inter_arrival = inter_arrival, type = type)
}
