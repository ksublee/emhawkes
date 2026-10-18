\name{exp_kernel_spec}
\alias{exp_kernel_spec}
\alias{as_exp_kernel}
\title{Numeric Exponential Kernel Specifications}
\description{Create a baseline-independent numeric kernel specification or extract
one from an existing specification. These interfaces are experimental.}
\usage{
exp_kernel_spec(alpha, beta, eta = NULL, type_col_map = NULL)
as_exp_kernel(object)
}
\arguments{
\item{alpha, beta}{Finite numeric matrices of equal positive dimensions, or
scalars representing one by one matrices. Decay rates must be positive and
equal within each row. Rows identify event types; columns identify kernels.}
\item{eta}{Optional finite numeric matrix of the same shape.}
\item{type_col_map}{One list entry per event type, containing distinct integer
kernel columns. Entries can overlap. When omitted, the number of columns must
equal the number of event types and event type i updates column i.}
\item{object}{An \code{hspec} with numeric kernel parameters and no impact
callback, or an existing \code{exp_kernel_spec}.}
}
\details{Construction and conversion execute no user functions and use no random
numbers. Function-valued kernel parameters and impact callbacks are rejected.
Baselines, mark generators, mark densities and residual laws are not copied.
Signed components and jumps are permitted; construction does not establish
that the corresponding model can be simulated. Execution state is separate
from the specification. Existing Hawkes calculation paths are unchanged.}
\value{An S3 \code{exp_kernel_spec} containing \code{alpha}, \code{beta},
\code{eta}, \code{dimens}, \code{kernel_cols} and the resolved
\code{type_col_map}.}
\seealso{\code{\link{exp_kernel_replay}}}
\examples{
h <- new("hspec", mu = .3, alpha = .2, beta = 1)
kernel <- as_exp_kernel(h)
exp_kernel_replay(kernel, c(0, .25, .5), lambda_component0 = .1)
}
