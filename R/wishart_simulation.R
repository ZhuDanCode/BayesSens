#' Simulate data from a 2-equation linear regression model with jointly Gaussian noise
#' @param n integer; number of datapoints.
#' @param p integer; number of measurements of first set of covariates.
#' @param k integer; number of measurements of second set of covariates.
#' @param beta (optional) numeric vector; the first set of regression coefficients.
#' It should be of length (p+1) if 'intercept_1 = T', and the intercept
#' should be the first element.
#' @param gamma (optional) numeric vector; the second set of regression coefficients.
#' It shoud be of length (k+2) if 'intercept_2 = T'. The intercept should be
#' the first element. The endogeneity factor should be the second element.
#' @param Sigma (optional) 2x2 numeric matrix; the covariance matrix of the distribution.
#' @param intercept_1 TRUE or FALSE; whether to include an intercept / bias term.
#' @param intercept_2 TRUE or FALSE; whether to include an intercept / bias term.
#' @return A list of the data X, y, Xs, s and the parameters beta, gamma and sigma.
#' @examples
#' data0 <- wishart_data(n = 100, p = 10, k = 3)
#' data1 <- wishart_data(n = 100, p = 2, k = 3,
#'   beta = c(0.2, -0.3, 0.4), gamma = c(0.2, 0.6, -0.2),
#'   Sigma = matrix(c(1, 0.2, 0.2, 1), 2, 2),
#'   intercept_1 = FALSE, intercept_2 = FALSE
#' )
#' @export
wishart_data <- function(n, p, k, beta, gamma, Sigma,
                         intercept_1 = TRUE, intercept_2 = TRUE) {

  num_covariates_1 <- ifelse(intercept_1, p + 2, p + 1)
  if (missing(beta)) beta <- rnorm(num_covariates_1)
  if (length(beta) != num_covariates_1)
    stop(paste("The length of beta does not match the number of covariates.",
               "Did you forget to specify the intercept or endogeneity?\n"))
  num_covariates_2 <- ifelse(intercept_2, k + 1, k)
  if (missing(gamma)) gamma <- rnorm(num_covariates_2)
  if (length(gamma) != num_covariates_2)
    stop(paste("The length of gamma does not match the number of covariates.",
               "Did you forget to specify the intercept?\n"))
  if (missing(Sigma)) Sigma <- pdmatrix(2)$Sigma

  err <- t(chol(Sigma)) %*% matrix(rnorm(2 * n), 2, n) %>% t()

  Xs <- matrix(rnorm(n * k), nrow = n, ncol = k)
  if (intercept_2) Xs <- cbind(1, Xs)
  s <- Xs %*% gamma + err[, 2]

  Xy <- matrix(rnorm(n * p), nrow = n, ncol = p)
  if (intercept_1) {
    Xy <- cbind(1, s, Xy)
  } else {
    Xy <- cbind(s, Xy)
  }
  y <- Xy %*% beta + err[, 1]

  list(Xy = Xy, y = y, beta = beta, Xs = Xs, s = s, gamma = gamma, Sigma = Sigma)
}
