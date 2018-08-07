#' #' Simulate data from a Vector AutoRegressive (VAR) model
#' #' @param m integer; number of datapoints.
#' #' @param n integer; dimension of the observation.
#' #' @param p integer; the time-series lag.
#' #' @param init_y_t numeric matrix (optional); the initial state of the time series.
#' #' The matrix must have as many rows as the lag, and the first
#' #' row should be the state at the first time point.
#' #' @param b_0 numeric vector (optional); the intercept.
#' #' @param B list of numeric matrices (optional); the regression coefficients.
#' #' @param Sigma numeric matrix (optional); the covariance matrix of the error.
#' #' @examples
#' #' simulate_VAR(30, 3, 2)
#' #' @export
#' simulate_VAR <- function(m, n, p, init_y_t, b_0, B, Sigma) {
#'   if (missing(b_0)) b_0 <- rnorm(n)
#'   if (missing(B)) B <- purrr::map(1:p, ~stationary_matrix(n, 1/p))
#'   if (missing(Sigma)) Sigma <- pdmatrix(n)$Sigma
#'
#'   y_t <- matrix(0, nrow = n, ncol = m)
#'   if (missing(init_y_t)) {
#'     y_t[,(m-p+1):m] <- matrix(rnorm(n*p), nrow = n, ncol = p)
#'   } else {
#'     assertthat::are_equal(ncol(init_y_t), p)
#'     y_t[,(m-p+1):m] <- t(init_y_t)
#'   }
#'
#'   for (i in (m-p):1) {
#'     y_t[,i] <- project_ts(b_0, B, y_t[, (i+1):(i+p), drop = FALSE], Sigma)
#'   }
#'
#'   list(data = t(y_t[, rev(1:m)]), model = list(b_0 = b_0, B = B, Sigma = Sigma))
#' }
#'
#' project_ts <- function(b, B, y_t, Sigma) {
#'   z <- b
#'   for (i in seq_along(B)) {
#'     z <- z + B[[i]] %*% y_t[,i]
#'   }
#'   z + MASS::mvrnorm(mu = numeric(length(b)), Sigma = Sigma)
#' }
#'
#' stationary_matrix <- function(k, upper_bound = 1) {
#'   res <- svd(pdmatrix(k)$Sigma)
#'   res$u %*% diag(upper_bound * res$d / max(res$d)) %*% t(res$v)
#' }
