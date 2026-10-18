library(emhawkes)
close <- function(x, y, tol = 1e-9) stopifnot(isTRUE(all.equal(x, y, tolerance = tol)))
close(ptzexp(.25, .5, 1), 29 / 132)
for (params in list(c(.5, 1), c(.5, .5), c(1, 1), c(3, 1))) {
  a <- params[1]; ell <- params[2]
  close(integrate(function(x) dtzexp(x, a, ell), 0, Inf,
                  subdivisions = 1000, rel.tol = 1e-10)$value, 1, tol = 1e-8)
  p <- c(0, 1e-8, .1, .5, .9, 1 - 1e-8)
  close(ptzexp(qtzexp(p, a, ell), a, ell), p)
  close(ptzexp(qtzexp(log(1 - p), a, ell, lower.tail = FALSE, log.p = TRUE), a, ell), p)
}
set.seed(428)
z <- rtzexp(50000, a = .5, ell = 1)
stopifnot(abs(mean(z <= .25) - ptzexp(.25, .5, 1)) < .01)
stopifnot(length(rtzexp(0, .5, 1)) == 0, length(rtzexp(1, .5, 1)) == 1)
stopifnot(inherits(tryCatch(dtzexp(4, 4, 1), error = identity), "error"))
stopifnot(is.finite(ptzexp(1000, .5, 1, lower.tail = FALSE, log.p = TRUE)))
close(qtzexp(ptzexp(1000, .5, 1, lower.tail = FALSE, log.p = TRUE),
             .5, 1, lower.tail = FALSE, log.p = TRUE), 1000)
tiny <- qtzexp(-1000, 1, sqrt(6), log.p = TRUE)
stopifnot(tiny > 0, is.finite(tiny))
close(ptzexp(tiny, 1, sqrt(6), log.p = TRUE), -1000)

# The ell=1, a->0 limit is unit exponential. Rescale the tiny lower-tail
# probabilities so an erroneous zero cannot pass an absolute-error check.
for (a in c(1e-12, 1e-16, 1e-100, 1e-300)) {
  close(ptzexp(qtzexp(c(.1, .5, .9), a, 1), a, 1), c(.1, .5, .9))
  close(ptzexp(a / 2, a, 1) / a, .5)
  close(qtzexp(ptzexp(a / 2, a, 1, log.p = TRUE), a, 1, log.p = TRUE) / a, .5)
}
