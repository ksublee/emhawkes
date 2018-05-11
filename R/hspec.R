setClassUnion("matrixORnumeric", c("matrix", "numeric"))
setClassUnion("functionOrNULL",members=c("function", "NULL"))

#' An S4 class to represent an exponential marked Hawkes model
#'
#' This class represents a specification of a marked Hawkes model with exponential kernel.
#' The intensity of the ground process is defined by:
#' \deqn{\lambda(t) = \mu + \int ( \alpha + \Psi ) * exp( -\beta (t-u)) d N(t)}
#' where \eqn{\alpha} represents exciting term that does not related to mark
#' \eqn{\Psi} represents the term related to mark and is represented by mark and impact.
#'
#' \code{mu}, \code{alpha} and \code{beta} are required slots for every exponential Hawkes model.
#' \code{mark} and \code{impact} are additional slots.
#'
#' @slot mu numeric value or matrix, automatically converted to matrix
#' @slot alpha numeric value or matrix, automatically converted to matrix, exciting term
#' @slot beta numeric value or matrix, automatically converted to matrix, exponential decay
#' @slot mark_hawkes a function that represets mark for counting process
#' @slot mark_lambda a function that represets mark for lambda process
#' @slot impact a function that describe the after impact of mark to lambda
#'
#' @examples
#' MU <- matrix(c(0.2), nrow = 2)
#' ALPHA <- matrix(c(0.75, 0.92, 0.92, 0.75), nrow = 2, byrow=TRUE)
#' BETA <- matrix(c(2.25, 2.25, 2.25, 2.25), nrow = 2, byrow=TRUE)
#' mhspec2 <- new("hspec", mu=MU, alpha=ALPHA, beta=BETA)
#'
#' @export
setClass("hspec",
  slots = list(
    mu = "matrixORnumeric",
    alpha = "matrixORnumeric",
    beta = "matrixORnumeric",
    mark_hawkes = "functionOrNULL",
    impact = "functionOrNULL"
  )
)

#' Initialize the hspec object
#'
#'
#' \code{mu}, \code{alpha} and \code{beta} are required slots for every exponential Hawkes model.
#' \code{mark} and \code{impact} are additional slots.
#'
#' @param .Object hspec
#' @param mu numeric value or matrix. Shoud be a matrix for two or larger dimensional model.
#' @param alpha numeric value or matrix. Shoud be a matrix for two or larger dimensional model.
#' @param beta numeric value or matrix. Shoud be a matrix for two or larger dimensional model.
#' @param mark_hawkes a function that generate mark
#' @param mark_lambda a function that represets mark for lambda process
#' @param impact a function that describe the mark impact
#' @param stability_check check the spectral radius
#'
#' @export
setMethod(
  "initialize",
  "hspec",
  function(.Object, mu, alpha, beta,
           mark_hawkes=NULL, mark_lambda=NULL, impact=NULL, stability_check=FALSE){

    # If mark_hawkes is not provided, then mark_hawkes is constant 1.
    if (is.null(mark_hawkes)) mark_hawkes <- function(...) 1

    # check the number of arguments
    if(length(formals(mark_hawkes)) == 0){
      .Object@mark_hawkes <- function(...) mark_hawkes()
    } else {
      .Object@mark_hawkes <- mark_hawkes
    }

    .Object@mu <- as.matrix(mu)
    .Object@alpha <- matrix(alpha, nrow = length(mu), ncol = length(mu))
    .Object@beta <- matrix(beta, nrow = length(mu), ncol = length(mu))
    .Object@impact <- impact

    # Check spectral radius, only works for non marked model
    if ( stability_check==TRUE && max(abs(eigen(alpha/beta)$values)) > 1)
      warning("This model may not satisfy the stability condition.")

    methods::callNextMethod()

    .Object


  }
)


