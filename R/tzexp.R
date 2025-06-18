#' Trapezoid + Exponential Distribution
#'
#' Density, distribution function, quantile function and random generation for a custom trapezoid + exponential distribution.
#'
#' @name tzexp
#' @aliases dtzexp ptzexp qtzexp rtzexp
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param a location parameter for transition (must be > 0)
#' @param ell rate parameter for exponential decay (must be > 0)
#' @return
#' \itemize{
#'   \item \code{dtzexp} gives the density
#'   \item \code{ptzexp} gives the distribution function
#'   \item \code{qtzexp} gives the quantile function
#'   \item \code{rtzexp} generates random deviates
#' }
#'
#' @description
#' These functions implement a custom distribution combining a trapezoidal section (0 < x < a)
#' and an exponential tail (x ≥ a). The distribution is parameterized by:
#' \itemize{
#'   \item \code{a}: transition point between trapezoid and exponential
#'   \item \code{ell}: rate parameter for the exponential tail
#' }
#'
#' @details
#' The trapezoid+exponential distribution has the probability density function:
#'
#' \deqn{
#' f(x) = \begin{cases}
#' 0 & \text{if } x \leq 0 \\
#' \frac{(p\ell - c)}{a} x + c & \text{if } 0 < x < a \\
#' p\ell e^{-\ell (x - a)} & \text{if } x \geq a
#' \end{cases}
#' }
#'
#' where:
#' \deqn{
#' p = \frac{\ell - \frac{a\ell}{3}}{\frac{a^2\ell^2}{6} + \frac{2a\ell}{3} + 1}
#' }
#' \deqn{
#' c = \frac{2 - 2p - p\ell a}{a}
#' }
#'
#' The trapezoid+exponential distribution has the following characteristics:
#' \itemize{
#'   \item Support on [0, ∞)
#'   \item Continuous probability distribution
#'   \item Linear density from 0 to a
#'   \item Exponential decay for x > a
#' }
#'
#'
#'
#' @examples
#' # Example usage of dtzexp function
#' x_values <- seq(-1, 5, by = 0.1)
#' a_value <- 2
#' ell_value <- 1
#' pdf_values <- dtzexp(x_values, a_value, ell_value)
#' plot(x_values, pdf_values, type = "l", main = "PDF of Custom Distribution")
#'
#' # Example usage of ptzexp function
#' cdf_values <- ptzexp(x_values, a_value, ell_value)
#' plot(x_values, cdf_values, type = "l", main = "CDF of Custom Distribution")
#'
#' # Example usage of rtzexp function
#' set.seed(123)
#' n_samples <- 1000
#' random_samples <- rtzexp(n_samples, a_value, ell_value)
#' hist(random_samples, breaks = 30, main = "Histogram of Random Samples")
#'
NULL

#' @rdname tzexp
#' @export
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

#' @rdname tzexp
#' @export
ptzexp <- function(q, a, ell){

  p <- (ell - a * ell/3) /(a^2 * ell^2 / 6 + 2 * a * ell / 3 + 1)

  c_ <- (2 - 2*p - p*ell*a) / a

  compute_cdf <- function(q){

    if (q <= 0){

      0

    } else if (q > 0 & q < a) {

      c_ * q - q^2 * (c_ - ell * p) / 2 / a


    } else {

      1 - p * exp(-ell * (q - a))

    }
  }

  sapply(q, compute_cdf)

}

#' @rdname tzexp
#' @export
qtzexp <- function(p, a, ell) {

  p_val <- (ell - a * ell / 3) / (a^2 * ell^2 / 6 + 2 * a * ell / 3 + 1)
  c_ <- (2 - 2 * p_val - p_val * ell * a) / a

  res <- numeric(length(p))

  res[p <= 0] <- 0
  res[p >= 1] <- Inf

  # 3. 0 < p <= (1 - p_val) i.e., 0 < x < a : quadratic
  idx_quad <- (p > 0) & (p <= 1 - p_val)
  if (any(idx_quad)) {
    p_quad <- p[idx_quad]

    A_val <- (c_ - ell * p_val) / (2 * a)
    discriminant <- c_^2 - 4 * A_val * p_quad
    discriminant <- pmax(discriminant, 0)  # non-negativity
    res[idx_quad] <- (2 * p_quad) / (c_ + sqrt(discriminant))
  }

  # 4. (1 - p_val) < p < 1 i.e., x >= a : exponential
  idx_exp <- (p > 1 - p_val) & (p < 1)
  if (any(idx_exp)) {
    p_exp <- p[idx_exp]
    res[idx_exp] <- a - (1 / ell) * log((1 - p_exp) / p_val)
  }

  res
}


#' @rdname tzexp
#' @export
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

