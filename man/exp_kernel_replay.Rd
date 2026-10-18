\name{exp_kernel_replay}
\alias{exp_kernel_replay}
\title{Replay an Observed Path Through a Numeric Exponential Kernel}
\description{Decay, integrate and update components along an observed path,
excluding the baseline, likelihood and random draws. This is an experimental
R reference interface, independent of the likelihood engine option.}
\usage{
exp_kernel_replay(kernel, inter_arrival, type = NULL, mark = NULL,
                  lambda_component0, history = FALSE)
}
\arguments{
\item{kernel}{An \code{exp_kernel_spec} or a convertible numeric \code{hspec}.}
\item{inter_arrival}{Finite nonnegative numeric vector beginning with zero.}
\item{type}{Event types, one per row. Required for multivariate kernels. Entries
after the initial row must be integer event types; the initial entry is ignored.}
\item{mark}{Finite numeric marks, one per row; defaults to zero initially and
one at each event. The initial mark is ignored.}
\item{lambda_component0}{Explicit finite initial intensity components excluding
the baseline. A scalar expands to every component. A matrix must match the kernel
shape; a vector of matching length uses R's column-major matrix order. Automatic
initialization policies belong to the calling model.}
\item{history}{Logical; retain full component histories?}
}
\details{The initial row causes no update. Zero-length event intervals apply
updates in input order. Negative components are not clipped. No terminal
observation interval is added. Count paths and initial counts are the calling
model's responsibility. Full component histories are allocated only when requested.}
\value{A list. For M rows including the initial row, D event types and K kernel
columns, \code{kernel_integral} and \code{kernel_before} are M by D matrices of
interval integrals and pre-event row sums. Initial integrals are zero; initial
sums are the initial component row sums. \code{final_component} is the D by K
state after the last event, or the initial state for a one-row path. With
\code{history = TRUE}, \code{component_before} and \code{component_after} are
M by (D times K), with row-major kernel ordering within each path row.}
\seealso{\code{\link{exp_kernel_spec}}}
\examples{
kernel <- exp_kernel_spec(.2, 1)
path <- exp_kernel_replay(kernel, c(0, .25), lambda_component0 = .1,
                          history = TRUE)
path$kernel_integral
path$final_component
}
