#' Numeric Exponential Kernel Calculations
#'
#' Create a kernel specification independent of a baseline or residual law, or
#' extract one from an existing numeric [hspec-class] object. These experimental
#' interfaces do not change the existing Hawkes calculation paths.
#'
#' @param alpha,beta Finite numeric matrices with the same positive dimensions,
#'   or numeric scalars (one by one matrices). `beta` must be positive and
#'   constant within each row. Rows identify event types; columns identify kernels.
#' @param eta Optional finite numeric matrix with the same shape as `alpha`.
#' @param type_col_map List of kernel columns to update for each event type.
#'   Columns must be distinct within each entry; overlap across entries is allowed.
#'   If omitted, the number of columns must equal the number of event types.
#' @param object An `hspec` with numeric `alpha`, `beta`, optional numeric `eta`
#'   and no `impact` callback, or an existing `exp_kernel_spec`.
#' @return A validated S3 `exp_kernel_spec` containing parameters, dimensions and
#'   the resolved column mapping. No execution state is stored in it.
#' @details Construction and conversion do not execute user functions or use
#'   random numbers. Function-valued kernel parameters and impact callbacks are
#'   rejected rather than evaluated or discarded. Baselines, mark generators,
#'   mark densities and residual laws belong to the calling model and are not
#'   copied. Signed components and jumps are permitted; this alone does not
#'   establish that a model can be simulated.
#' @seealso [exp_kernel_replay()]
#' @examples
#' h <- new("hspec", mu = .3, alpha = .2, beta = 1)
#' kernel <- as_exp_kernel(h)
#' exp_kernel_replay(kernel, c(0, .25, .5), lambda_component0 = .1)
#' @export
exp_kernel_spec <- function(alpha, beta, eta = NULL, type_col_map = NULL) {
  numeric_matrix <- function(value, name) {
    if (!is.numeric(value) || is.object(value) || !length(value) ||
        any(!is.finite(value)) ||
        (!is.null(dim(value)) && length(dim(value)) != 2L) ||
        (is.null(dim(value)) && length(value) != 1L))
      stop(name, " must be a finite numeric matrix or scalar.", call. = FALSE)
    if (is.null(dim(value))) value <- matrix(value, 1L, 1L)
    value
  }
  alpha <- numeric_matrix(alpha, "alpha")
  beta <- numeric_matrix(beta, "beta")
  if (!identical(dim(alpha), dim(beta)))
    stop("alpha and beta must have the same dimensions.", call. = FALSE)
  if (any(beta <= 0)) stop("beta must be positive.", call. = FALSE)
  if (any(beta != beta[, 1L]))
    stop("The numeric kernel interface requires equal beta values within each row.", call. = FALSE)
  if (!is.null(eta)) {
    eta <- numeric_matrix(eta, "eta")
    if (!identical(dim(eta), dim(alpha)))
      stop("eta must match the kernel matrix dimensions.", call. = FALSE)
  }
  dimens <- nrow(alpha)
  kernel_cols <- ncol(alpha)
  if (is.null(type_col_map)) {
    if (kernel_cols != dimens)
      stop("type_col_map is required when kernel columns differ from event types.", call. = FALSE)
    type_col_map <- as.list(seq_len(dimens))
  }
  if (!is.list(type_col_map) || length(type_col_map) != dimens ||
      !all(vapply(type_col_map, function(columns) {
        is.numeric(columns) && !is.object(columns) && is.null(dim(columns)) &&
          length(columns) > 0L && all(is.finite(columns)) &&
          all(columns == floor(columns) & columns >= 1 & columns <= kernel_cols) &&
          !anyDuplicated(columns)
      }, logical(1))))
    stop("type_col_map must contain one nonempty list entry per event type, with distinct integer kernel columns in range.", call. = FALSE)
  structure(list(alpha = alpha, beta = beta, eta = eta,
    dimens = dimens, kernel_cols = kernel_cols,
    type_col_map = lapply(type_col_map, as.integer)), class = "exp_kernel_spec")
}

#' @rdname exp_kernel_spec
#' @export
as_exp_kernel <- function(object) {
  if (inherits(object, "exp_kernel_spec"))
    return(exp_kernel_spec(object$alpha, object$beta, object$eta, object$type_col_map))
  if (!methods::is(object, "hspec"))
    stop("object must be an hspec or exp_kernel_spec.", call. = FALSE)
  if (!is.null(object@impact))
    stop("The numeric kernel interface does not support impact callbacks.", call. = FALSE)
  # Read slots directly: setting() can execute parameter callbacks.
  exp_kernel_spec(object@alpha, object@beta, object@eta, object@type_col_map)
}

# These helpers take already validated inputs. Their state is local to one run.
exp_kernel_advance <- function(kernel, component, dt) {
  # The constructor guarantees equal rates within each row. R recycles these
  # row weights down each column, avoiding duplicate decay/integral work.
  rates <- kernel$beta[, 1L]
  list(before = component * exp(-rates * dt),
       integral = integrate_exp_kernel(component, rates, dt))
}

exp_kernel_update <- function(kernel, component, event_type, mark) {
  columns <- kernel$type_col_map[[event_type]]
  # Retain the existing order: alpha addition, then marked eta addition. Keep
  # additions on unselected columns too (including the signed-zero behavior).
  jump <- matrix(0, kernel$dimens, kernel$kernel_cols)
  jump[, columns] <- kernel$alpha[, columns]
  component <- component + jump
  if (!is.null(kernel$eta)) {
    jump[,] <- 0
    jump[, columns] <- kernel$eta[, columns] * (mark - 1)
    component <- component + jump
  }
  component
}

#' Replay an Observed Path Through a Numeric Exponential Kernel
#'
#' Decay, integrate and update kernel components along an observed path. This
#' deterministic calculation excludes the baseline, likelihood and random draws.
#'
#' @param kernel An [exp_kernel_spec()] or a convertible numeric `hspec`.
#' @param inter_arrival Finite nonnegative numeric vector beginning with zero.
#' @param type Event types, one per row. Required for multivariate kernels.
#'   After the initial row, entries must be integers from 1 to the number of rows
#'   in `alpha`. The initial entry is ignored.
#' @param mark Finite numeric marks, one per row; defaults to zero initially and
#'   one at each event. The initial mark is ignored.
#' @param lambda_component0 Explicit finite initial intensity components, excluding
#'   the baseline. A scalar expands to every component; a matrix must have the
#'   kernel shape. A vector of matching length uses R's column-major matrix order.
#'   Automatic initialization policies are handled by the calling model.
#' @param history Logical; retain full pre-event and post-event component histories?
#' @return A list with `kernel_integral` and `kernel_before`, each M by D, and
#'   `final_component`, D by K. Here M includes the initial row, D is the number
#'   of event types and K the kernel columns. Initial integrals are zero and
#'   initial pre-event sums equal the initial component row sums. If `history`
#'   is true, `component_before` and `component_after` are M by (D times K), with
#'   each path row flattened in row-major kernel order. No baseline is included.
#' @details The initial row causes no event update. Zero-length event intervals
#'   still apply updates in input order. Negative components are not clipped.
#'   The final state is the state after the last recorded event; no terminal
#'   observation interval is added. Count paths and initial counts remain the
#'   responsibility of the calling model. With `history = FALSE`, full component
#'   histories are not allocated. This experimental R reference interface is
#'   independent of the `emhawkes.engine` likelihood option.
#' @seealso [exp_kernel_spec()]
#' @examples
#' kernel <- exp_kernel_spec(.2, 1)
#' path <- exp_kernel_replay(kernel, c(0, .25), lambda_component0 = .1,
#'                           history = TRUE)
#' path$kernel_integral
#' path$final_component
#' @export
exp_kernel_replay <- function(kernel, inter_arrival, type = NULL, mark = NULL,
                              lambda_component0, history = FALSE) {
  kernel <- as_exp_kernel(kernel)
  if (!is.logical(history) || length(history) != 1L || is.na(history))
    stop("history must be TRUE or FALSE.", call. = FALSE)
  size <- length(inter_arrival)
  if (!is.numeric(inter_arrival) || is.object(inter_arrival) ||
      !is.null(dim(inter_arrival)) || !size || any(!is.finite(inter_arrival)) ||
      any(inter_arrival < 0) || inter_arrival[1L] != 0)
    stop("inter_arrival must be a finite nonnegative vector starting with zero.", call. = FALSE)
  if (is.null(type)) {
    if (kernel$dimens != 1L)
      stop("type is required for a multivariate kernel.", call. = FALSE)
    type <- c(0, rep(1, size - 1L))
  }
  if (!is.numeric(type) || is.object(type) || !is.null(dim(type)) ||
      length(type) != size || any(!is.finite(type)) ||
      any(type[-1L] != floor(type[-1L]) | type[-1L] < 1 | type[-1L] > kernel$dimens))
    stop("type must match inter_arrival and contain valid integer event types after the initial row.", call. = FALSE)
  if (is.null(mark)) mark <- c(0, rep(1, size - 1L))
  if (!is.numeric(mark) || is.object(mark) || !is.null(dim(mark)) ||
      length(mark) != size || any(!is.finite(mark)))
    stop("mark must be finite and match inter_arrival.", call. = FALSE)
  if (missing(lambda_component0) || is.null(lambda_component0) ||
      !is.numeric(lambda_component0) || is.object(lambda_component0))
    stop("Supply explicit finite numeric initial intensity components in lambda_component0.", call. = FALSE)
  component <- normalize_lambda_initial(lambda_component0, kernel$dimens, kernel$kernel_cols)
  integral <- before <- matrix(0, size, kernel$dimens)
  before[1L, ] <- rowSums(component)
  if (history) {
    component_before <- component_after <- matrix(0, size, length(component))
    component_before[1L, ] <- component_after[1L, ] <- as.vector(t(component))
  }
  for (n in seq_len(size - 1L) + 1L) {
    advanced <- exp_kernel_advance(kernel, component, inter_arrival[n])
    integral[n, ] <- rowSums(advanced$integral)
    before[n, ] <- rowSums(advanced$before)
    component <- exp_kernel_update(kernel, advanced$before, type[n], mark[n])
    if (history) {
      component_before[n, ] <- as.vector(t(advanced$before))
      component_after[n, ] <- as.vector(t(component))
    }
  }
  result <- list(kernel_integral = integral, kernel_before = before,
                 final_component = component)
  if (history) {
    result$component_before <- component_before
    result$component_after <- component_after
  }
  result
}
