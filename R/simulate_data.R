#' Simulate data from a linear regression model with gaussian noise
#' @param n integer; number of datapoints.
#' @param p integer; number of measurements.
#' @param beta (optional) numeric vector; the regression coefficients.
#' @param sigma positive number; the standard deviation of the distribution.
#' @return A list of the data X, y, and the parameters beta, sigma.
#' @examples
#' data0 <- gaussian_data(n = 100, p = 10)
#' data1 <- gaussian_data(n = 100, p = 3, beta = c(0.2, -0.3, 0.4), sigma = 2)
#' @export
gaussian_data <- function(n, p, beta, sigma) {
  if (missing(beta)) beta <- rnorm(p)
  if (missing(sigma)) sigma <- runif(1)

  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  y <- X %*% beta + rnorm(n, sd = sigma)
  list(X = X, y = y, beta = beta, sigma = sigma)
}
