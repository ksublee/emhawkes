#' @include hspec.R
NULL

setGeneric("get_lambda0", function(object, ...) standardGeneric("get_lambda0"))
# Get the long-run expectation of intensities - lambda matrix
#
# This function is not crucial but useful when you want to set lamba0 by default.
#
#' @param object hspec
setMethod(
  f = "get_lambda0",
  signature(object = "hspec"),
  definition = function(object, inter_arrival = inter_arrival, type = type, mark = mark,
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
        mu0 <- mu(n = 1, inter_arrival = inter_arrival, type = type, mark = mark,
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


#' Compute Hawkes volatility
#'
#' This function computes Hawkes volatilty.
#'
#' @param object hspec
#' @param horizon
#' @param inter_arrival
#' @param type
#' @param mark
#' @param dependece
#' @param lambda0
#'
#' @rdname hvol
#' @export
setGeneric("hvol", function(object, horizon = 1,
                           inter_arrival = NULL,
                           type = NULL,
                           mark = NULL,
                           dependence = FALSE,
                           lambda0 = NULL, ...) standardGeneric("hvol"))
#'
#'
#' @rdname hvol
setMethod(
  f = "hvol",
  signature(object = "hspec"),
  definition = function(object, horizon = 1,
                        inter_arrival = NULL,
                        type = NULL,
                        mark = NULL,
                        dependence = FALSE,
                        lambda0 = NULL,
                        ...){

    h <- object

    infer_result <- infer_lambda(h, inter_arrival = inter_arrival, type = type, mark = mark, lambda0 = lambda0)

    Mu <- matrix(h@mu, nrow=2)
    Alpha <- h@alpha
    Beta <- diag(diag(h@beta))
    Eta <- h@eta

    up_move_idx <- type == 1
    down_move_idx <- type == 2


    #### compute Hawkes vol

    if (dependence){

      Z_1 <- sum(mark[up_move_idx] * infer_result$lambda[up_move_idx,1]) /
        sum(infer_result$lambda[up_move_idx,1])

      Z_2 <- sum(mark[down_move_idx] * infer_result$lambda[down_move_idx,2]) /
        sum(infer_result$lambda[down_move_idx,2])

      Z_sq_1 <- sum(mark[up_move_idx]^2 * infer_result$lambda[up_move_idx,1]) /
        sum(infer_result$lambda[up_move_idx,1])

      Z_sq_2 <- sum(mark[down_move_idx]^2 * infer_result$lambda[down_move_idx,2]) /
        sum(infer_result$lambda[down_move_idx,2])

      Z_ln_11 <- sum(mark[up_move_idx] * infer_result$lambda[up_move_idx,1] * infer_result$Nc[up_move_idx,1]) /
        sum(infer_result$lambda[up_move_idx,1] * infer_result$Nc[up_move_idx,1])

      Z_ln_12 <- sum(mark[down_move_idx] * infer_result$lambda[down_move_idx,1] * infer_result$Nc[down_move_idx,2]) /
        sum(infer_result$lambda[down_move_idx,1] * infer_result$Nc[down_move_idx,2])

      Z_ln_21 <- sum(mark[up_move_idx] * infer_result$lambda[up_move_idx,2] * infer_result$Nc[up_move_idx,1]) /
        sum(infer_result$lambda[up_move_idx,2] * infer_result$Nc[up_move_idx,1])

      Z_ln_22 <- sum(mark[down_move_idx] * infer_result$lambda[down_move_idx,2] * infer_result$Nc[down_move_idx,2]) /
        sum(infer_result$lambda[down_move_idx,2] * infer_result$Nc[down_move_idx,2])


      Z_nl_11 <- sum(mark[up_move_idx] * infer_result$Nc[up_move_idx,1] * infer_result$lambda[up_move_idx,1]) /
        sum(infer_result$Nc[up_move_idx,1] * infer_result$lambda[up_move_idx,1])

      Z_nl_12 <- sum(mark[down_move_idx] * infer_result$Nc[down_move_idx,1] * infer_result$lambda[down_move_idx,2]) /
        sum(infer_result$Nc[down_move_idx,1] * infer_result$lambda[down_move_idx,2])

      Z_nl_21 <- sum(mark[up_move_idx] * infer_result$Nc[up_move_idx,2] * infer_result$lambda[up_move_idx,1]) /
        sum(infer_result$Nc[up_move_idx,2] * infer_result$lambda[up_move_idx,1])

      Z_nl_22 <- sum(mark[down_move_idx] * infer_result$Nc[down_move_idx,2] * infer_result$lambda[down_move_idx,2]) /
        sum(infer_result$Nc[down_move_idx,2] * infer_result$lambda[down_move_idx,2])


      Z_ll_11 <- sum(mark[up_move_idx] * infer_result$lambda[up_move_idx,1] * infer_result$lambda[up_move_idx,1]) /
        sum(infer_result$lambda[up_move_idx,1] * infer_result$lambda[up_move_idx,1])

      Z_ll_12 <- sum(mark[down_move_idx] * infer_result$lambda[down_move_idx,1] * infer_result$lambda[down_move_idx,2]) /
        sum(infer_result$lambda[down_move_idx,1] * infer_result$lambda[down_move_idx,2])

      Z_ll_21 <- sum(mark[up_move_idx] * infer_result$lambda[up_move_idx,2] * infer_result$lambda[up_move_idx,1]) /
        sum(infer_result$lambda[up_move_idx,2] * infer_result$lambda[up_move_idx,1])

      Z_ll_22 <- sum(mark[down_move_idx] * infer_result$lambda[down_move_idx,2] * infer_result$lambda[down_move_idx,2]) /
        sum(infer_result$lambda[down_move_idx,2] * infer_result$lambda[down_move_idx,2])


      Z <- matrix(c(Z_1, Z_2,
                    Z_1, Z_2), nrow=2, byrow=TRUE)

      Z_sq <- matrix(c(Z_sq_1, Z_sq_2,
                       Z_sq_1, Z_sq_2), nrow=2, byrow=TRUE)

      Z_ln <- matrix(c(Z_ln_11, Z_ln_12,
                       Z_ln_21, Z_ln_22), nrow=2, byrow=TRUE)

      Z_nl <- matrix(c(Z_nl_11, Z_nl_12,
                       Z_nl_21, Z_nl_22), nrow=2, byrow=TRUE)

      Z_ll <- matrix(c(Z_ll_11, Z_ll_12,
                       Z_ll_21, Z_ll_22), nrow=2, byrow=TRUE)

      Z_nl <- Z


      E_lambda <- solve(Beta - Alpha + Eta - Eta * Z) %*% Beta %*% Mu

      # To find E_ll : E_lambda_lambda

      G <- (Alpha - Eta + Eta * Z) %*% diag(c(E_lambda)) %*% t(Alpha - Eta) +
        (Alpha - Eta) %*% diag(c(E_lambda)) %*% t(Eta * Z) +
        (Eta * Z_sq^(0.5)) %*% diag(c(E_lambda)) %*% t(Eta * Z_sq^(0.5))


      H <- kronecker(diag(2), Alpha - Beta) + kronecker(Alpha - Beta, diag(2)) +
        rbind(cbind(Eta %*% diag(t(Z_ll)[,1] - 1), matrix(rep(0,4), nrow=2)),
              cbind(matrix(rep(0,4), nrow=2), Eta %*% diag(t(Z_ll)[,2] - 1))) +
        rbind(cbind(t(Eta)[,1] %*% diag(t(Z_ll)[,1] - 1), matrix(rep(0,2), nrow=1)),
              cbind(matrix(rep(0,2), nrow=1), t(Eta)[,1] %*% diag(t(Z_ll)[,2] - 1)),
              cbind(t(Eta)[,2] %*% diag(t(Z_ll)[,1] - 1), matrix(rep(0,2), nrow=1)),
              cbind(matrix(rep(0,2), nrow=1), t(Eta)[,2] %*% diag(t(Z_ll)[,2] - 1)))


      C <- Beta %*% Mu %*% t(E_lambda) + t(Beta %*% Mu  %*% t(E_lambda)) + G

      E_ll <- matrix(solve(H, -c(C)), nrow=2)

      HA <- kronecker(diag(2), Alpha - Beta) +
        rbind(cbind(Eta %*% diag(t(Z)[,1] - 1), matrix(rep(0,4), nrow=2)),
              cbind(matrix(rep(0,4), nrow=2), Eta %*% diag(t(Z)[,2] - 1)))

      CA <- (Beta %*% Mu) %*% t(E_lambda) %*% diag(diag(Z))

      A <- t(matrix(solve(HA, -c(CA)), nrow=2))
      # A <- diag(diag(Z)) %*% E_lambda %*% t(E_lambda)

      CB <- Z_ll * E_ll + ((Alpha - Eta) * Z + Eta * Z_sq) %*% diag(c(E_lambda)) - t(A)

      B <- t(matrix(solve(HA, -c(CB)), nrow=2))

      E_nl <- A * horizon + B

      P1 <- Z_nl * 0.5 * A + t(Z_nl * 0.5 * A)
      P2 <- Z_nl * B  + t(Z_nl * B)
      P3 <- Z_sq * diag(c(E_lambda))

      E_NN <- P1 * horizon^2 + P2 * horizon + P3 * horizon

      variance <- matrix(c(1, -1), nrow=1) %*% (P2 * horizon + P3 * horizon) %*% matrix(c(1, -1))

      (std <- sqrt(variance))

    } else {

      ##### independent

      Z_ind_1 <- mean(mark[up_move_idx])
      Z_ind_2 <- mean(mark[down_move_idx])

      Z_sq_ind_1 <- mean(mark[up_move_idx]^2)
      Z_sq_ind_2 <- mean(mark[down_move_idx]^2)

      Z_ind <- matrix(c(Z_ind_1, Z_ind_2,
                        Z_ind_1, Z_ind_2), nrow=2, byrow=TRUE)

      Z_sq_ind <-  matrix(c(Z_sq_ind_1, Z_sq_ind_2,
                            Z_sq_ind_1, Z_sq_ind_2), nrow=2, byrow=TRUE)

      E_lambda_ind <- solve(Beta - Alpha + Eta - Eta * Z_ind) %*% Beta %*% Mu

      G_ind <- (Alpha - Eta + Eta * Z_ind) %*% diag(c(E_lambda_ind)) %*% t(Alpha - Eta) +
        (Alpha - Eta) %*% diag(c(E_lambda_ind)) %*% t(Eta * Z_ind) +
        (Eta * Z_sq_ind^(0.5)) %*% diag(c(E_lambda_ind)) %*% t(Eta * Z_sq_ind^(0.5))

      H_ind <- kronecker(diag(2), Alpha - Beta) + kronecker(Alpha - Beta, diag(2)) +
        rbind(cbind(Eta %*% diag(t(Z_ind)[,1] - 1), matrix(rep(0,4), nrow=2)),
              cbind(matrix(rep(0,4), nrow=2), Eta %*% diag(t(Z_ind)[,2] - 1))) +
        rbind(cbind(t(Eta)[,1] %*% diag(t(Z_ind)[,1] - 1), matrix(rep(0,2), nrow=1)),
              cbind(matrix(rep(0,2), nrow=1), t(Eta)[,1] %*% diag(t(Z_ind)[,2] - 1)),
              cbind(t(Eta)[,2] %*% diag(t(Z_ind)[,1] - 1), matrix(rep(0,2), nrow=1)),
              cbind(matrix(rep(0,2), nrow=1), t(Eta)[,2] %*% diag(t(Z_ind)[,2] - 1)))

      C_ind <- Beta %*% Mu %*% t(E_lambda_ind) + t(Beta %*% Mu  %*% t(E_lambda_ind)) + G_ind

      E_ll_ind <- matrix(solve(H_ind, -c(C_ind)), nrow=2)

      HA_ind <- kronecker(diag(2), Alpha - Beta) +
        rbind(cbind(Eta %*% diag(t(Z_ind)[,1] - 1), matrix(rep(0,4), nrow=2)),
              cbind(matrix(rep(0,4), nrow=2), Eta %*% diag(t(Z_ind)[,2] - 1)))

      CA_ind <- (Beta %*% Mu) %*% t(E_lambda_ind) %*% diag(diag(Z_ind))

      A_ind <- t(matrix(solve(HA_ind, -c(CA_ind)), nrow=2))

      CB_ind <- Z_ind * E_ll_ind + ((Alpha - Eta) * Z_ind + Eta * Z_sq_ind) %*% diag(c(E_lambda_ind)) - t(A_ind)

      B_ind <- t(matrix(solve(HA_ind, -c(CB_ind)), nrow=2))

      E_nl_ind <- A_ind * horizon + B_ind

      P1_ind <- Z_ind * 0.5 * A_ind + t(Z_ind * 0.5 * A_ind)
      P2_ind <- Z_ind * B_ind  + t(Z_ind * B_ind)
      P3_ind <- diag(diag(Z_sq_ind)) %*% diag(c(E_lambda_ind))

      E_NN_ind <- P1_ind * horizon^2 + P2_ind * horizon + P3_ind * horizon

      variance_ind <- matrix(c(1, -1), nrow=1) %*% (P2_ind * horizon + P3_ind * horizon) %*% matrix(c(1, -1))

      (std_ind <- sqrt(variance_ind))

    }

  }
)

