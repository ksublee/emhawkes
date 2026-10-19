\name{exp_kernel_plan}
\alias{exp_kernel_plan}
\title{Prepare Reusable Numeric Exponential Kernel Calculations}
\usage{
exp_kernel_plan(kernel, mu)
}
\arguments{
\item{kernel}{An exp_kernel_spec or convertible numeric hspec.}
\item{mu}{Finite nonnegative baseline, scalar or one value per event type.}
}
\value{
A locked environment with two functions:
\code{step(component, dt, type = NULL, mark = 1, integrate = TRUE)} and
\code{wait(component, residual)}. Their results match
\code{\link{exp_kernel_step}} and \code{\link{exp_kernel_wait}}.
}
\description{
Validate and capture a numeric exponential kernel and baseline once for
repeated deterministic calculations.
}
\details{
Later edits to the original inputs do not change the plan. Every call still
checks the current component and changing arguments. Components are passed
explicitly; no current state, histories or random generator are stored or
advanced. Create a new plan when kernel parameters or baseline change.
This experimental numeric interface requires equal beta values within each
row. The likelihood engine option and existing event loops are unchanged.
}
\examples{
kernel <- exp_kernel_spec(.2, 1)
plan <- exp_kernel_plan(kernel, mu = .3)
component <- matrix(.1)
dt <- plan$wait(component, residual = .7)
next_state <- plan$step(component, dt, type = 1)
next_state$after
}
