\name{parameter_plan}
\alias{parameter_plan}
\title{Register and Restore Model Parameters}
\usage{
parameter_plan(components, prefixes = names(components),
               named = character(), reduced = TRUE)
}
\arguments{
\item{components}{Named list of numeric matrices/scalars, functions or NULL.}
\item{prefixes}{Matrix name prefixes, a character vector matching components.}
\item{named}{Names of explicit named parameter-vector components. Their names
are retained; equal starting values do not imply sharing.}
\item{reduced}{Share equal starting values within each matrix?}
}
\value{
A locked environment with \code{start}, the full starting vector, and
\code{restore(param)}, returning a named list of matrices, functions and vectors.
}
\description{
Prepare model-independent parameter coordinates and component reconstruction.
}
\details{
Components are visited in list order. Matrix naming and sharing match
\code{\link{hfit}}. Function parameters come from their named numeric
\code{param} defaults, evaluated in the original environment; function bodies
are not called. Names shared across components identify one parameter.
Conflicting starts warn and use the first value. Functions retain their
environments and binding semantics. No model, optimizer, initialization policy
or random generator is supplied by this API. External function environments
remain references; evaluating user parameter defaults can have side effects.

Restore accepts the full vector in starting order, or exactly its names in
any order. Nonfinite proposals are retained for the objective to handle.
Fixed/free selectors and constraints use this full coordinate order and are
handled by the optimizer layer. Create a new plan when definitions change.
}
\examples{
plan <- parameter_plan(list(mu = matrix(c(.2, .2), 2),
  residual = c(shape = 2, rate = 2)), named = "residual")
plan$start
plan$restore(plan$start)
}
