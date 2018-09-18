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


#' Frobenius norm of a matrix
#' @param A A numeric matrix
#' @export
frobenius_norm <- function(A) { sqrt(sum(A^2)) }


#' Maximum norm of a matrix
#' @param A A numeric matrix
#' @export
maximum_norm <- function(A) { max(abs(A)) }


#' Remove the burn-in samples
#' @param res The Jacobian array.
#' @param num_burns Positive integer; the number of samples to burn.
#' @export
burn_ins <- function(res, num_burns) {
  res[-seq(num_burns), ]
}


#' Thinning the samples
#' @param vec0 A numeric vector; the samples.
#' @param every_n Positive integer; take one sample every_n sample.
thinning <- function(vec0, every_n) {
  vec0[seq(1, length(vec0), every_n)]
}
