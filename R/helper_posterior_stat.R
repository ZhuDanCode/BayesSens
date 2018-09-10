#' Cumulative summary statistics
#' @param x numeric vector.
#' @param fun function; the summary statistics function.
#' @param ... Extra parameters to be passed to stat_fun.
#' @export
cumulative_stat <- function(x, fun, ...) {
  end_pt <- seq_along(x)
  f <- function(y) { fun(y, ...) }
  purrr::map_dbl(end_pt, ~f(x[1:.x]))
}


#' Jacobian matrix
#' @param x the sensitivity of the posterior samples
#' @param stat_fun function; the summary statistics function,
#' default to be the posterior mean.
#' @param ... extra parameters to be passed to `stat_fun`.
#' @export
jacobian <- function(x, stat_fun = mean, ...) {
  f <- function(y) { stat_fun(y, ...) }
  apply(x, c(2, 3), f)
}
