## probability functions for example

dtzexp <- function(x, a, ell){


  p <- (ell - a * ell/3) /(a^2 * ell^2 / 6 + 2 * a * ell / 3 + 1)

  c_ <- (2 - 2*p - p*ell*a) / a

  compute_pdf <- function(x){

    if (x > 0 & x < a) {

      (p*ell - c_) / a * x + c_

    } else if (x <= 0 ){

      0

    } else {

      p * ell * exp(-ell * (x - a))

    }
  }

  sapply(x, compute_pdf)

}

ptzexp <- function(x, a, ell){

  p <- (ell - a * ell/3) /(a^2 * ell^2 / 6 + 2 * a * ell / 3 + 1)

  c_ <- (2 - 2*p - p*ell*a) / a

  compute_cdf <- function(x){

    if (x <= 0){

      0

    } else if (x > 0 & x < a) {

      c_ * x - x^2 * (c_ - ell * p) / 2 / a


    } else {

      1 - p * exp(-ell * (x - a))

    }
  }

  sapply(x, compute_cdf)

}

rtzexp <- function(n, a, ell){


  p <- (ell - a * ell/3) /(a^2 * ell^2 / 6 + 2 * a * ell / 3 + 1)

  c_ <- (2 - 2*p - p*ell*a) / a


  ## first selection

  idx_exp <- rbinom(n, size = 1, prob = p)

  rv_exp <- rexp(sum(idx_exp == 1), ell) + a

  n_to_generate_trapziod <- sum(idx_exp == 0)

  while(TRUE){

    # 3 times are enough?
    X <- runif(n_to_generate_trapziod * 3, 0, a)
    Y <- runif(n_to_generate_trapziod * 3, 0, p*ell)

    selected_X <- X[Y < dtzexp(X, a, ell)]
    if(length(selected_X) >= n_to_generate_trapziod) break

  }
  rv_trapziod <- selected_X[1:n_to_generate_trapziod]

  result <- rep(0, n)
  result[which(idx_exp == 1)] <- rv_exp
  result[which(idx_exp == 0)] <- rv_trapziod

  result

}
