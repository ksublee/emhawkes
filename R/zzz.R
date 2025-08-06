# R/zzz.R

# Declare global variables to avoid NOTE about no visible binding during R CMD check
utils::globalVariables(c(
  "mu0", "alpha0", "beta0", "eta0", "impact0",
  "dmark0", "dresidual0", "presidual0", "qresidual0",
  "mu", "alpha", "beta", "eta"
))

.onLoad <- function(libname, pkgname) {
  # Code to run when the package is loaded can be added here
  invisible()
}
