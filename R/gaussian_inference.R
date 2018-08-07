#' Model inference for normal regression model with multivariate normal
#' prior for the mean and inverse gamma for the variance.
#' @param X A numeric matrix; the covariates.
#' @param y A numeric vector; the response variable.
#' @param b_0 A numeric vector; the mean for the multivariate normal prior.
#' @param B_0 A numeric matrix; the covariance for the multivariate normal prior.
#' @param alpha_0 A positive number; the shape parameter of the inverse-gamma prior.
#' @param delta_0 A positive number; the rate parameter of the inverse-gamma prior.
#' @param init_sigma (Optional) A positive number; the starting value of the sigma parameter.
#' @param num_steps Integer; number of MCMC steps.
#' @return A list of posterior samples of the parameters; each parameter is either a
#' a vector of length `num_steps`, or a matrix of dimension `num_steps` x d, where
#' d is the dimension of the parameter.
#' @examples
#' \dontrun{
#' n <- 100
#' p <- 10
#' data0 <- gaussian_data(n, p, intercept = T)
#' res <- gaussian_Gibbs(
#'   X = data0$X, y = data0$y,
#'   b_0 = rnorm(p+1), B_0 = pdmatrix(p+1)$Sigma,  # add one for the intercept
#'   alpha_0 = 13, delta_0 = 8
#' )
#' }
#' @export
gaussian_Gibbs <- function(X, y, b_0, B_0, alpha_0, delta_0,
                           init_sigma, num_steps = 1e4) {
  if (missing(init_sigma))
    init_sigma <- 1 / sqrt(rgamma(1, alpha_0 / 2, delta_0 / 2))

  # Initialisation
  n <- length(y)
  alpha_1 <- alpha_0 + n
  sigma_g <- init_sigma
  inv_B_0 <- solve(B_0)
  inv_B_0_times_b_0 <- inv_B_0 %*% b_0
  XTX <- crossprod(X)
  XTy <- crossprod(X, y)
  beta_res <- matrix(0, num_steps, length(b_0))
  sigma_res <- numeric(num_steps)

  pb <- txtProgressBar(1, num_steps, style = 3)
  for (i in 1:num_steps) {
    # Update beta
    B_g <- solve(sigma_g^(-2) * XTX + inv_B_0)
    b_g <- B_g %*% (sigma_g^(-2) * XTy + inv_B_0_times_b_0)
    beta_g <- MASS::mvrnorm(1, b_g, B_g)

    # Update sigma
    delta_g <- delta_0 + sum((y - X %*% beta_g)^2)
    sigma_g <- 1 / sqrt(rgamma(1, alpha_1 / 2, delta_g / 2))

    # Keep track
    beta_res[i,] <- beta_g
    sigma_res[i] <- sigma_g
    setTxtProgressBar(pb, i)
  }

  list(sigma = sigma_res, beta = beta_res)
}
