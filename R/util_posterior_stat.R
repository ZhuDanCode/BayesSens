#' Cumulative summary statistics
#' @param x numeric vector.
#' @param stat_fun function; the summary statistics function.
#' @param ... Extra parameters to be passed to stat_fun.
#' @export
cumulative_stat <- function(x, stat_fun, ...) {
  end_pt <- seq_along(x)
  purrr::map_dbl(end_pt, ~stat_fun(x[1:.x], ...))
}


#' Jacobian matrix
#' @param x the sensitivity of the posterior samples
#' @param dim0 the dimension of the Jacobian matrix
#' @export
jacobian <- function(x, dim0) {
  res <- colMeans(x)
  dim(res) <- dim0
  res
}
