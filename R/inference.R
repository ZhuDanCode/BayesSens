#' Model inference for normal regression model with multivariate normal
#' prior for the mean and inverse gamma for the variance.
#' @param X A numeric matrix; the covariates.
#' @param y A numeric vector; the response variable.
#' @param b_0 A numeric vector; the mean for the multivariate normal prior.
#' @param B_0 A numeric matrix; the covariance for the multivariate normal prior.
#' @param alpha_0 A positive number; the shape parameter of the inverse-gamma prior.
#' @param delta_0 A positive number; the rate parameter of the inverse-gamma prior.
#' @param beta (Optional) A numeric vector
#' @param sigma (Optional) A numeric matrix
#' @param num_steps Integer; number of MCMC steps.
#' @param burn_ins Integer; number of burn-ins.
#' @examples
#' \dontrun{
#' n <- 100
#' p <- 10
#' data0 <- gaussian_data(n, p)
#' res <- gaussian_Gibbs(data0$X, data0$y,
#'   b_0 = rnorm(p), B_0 = pdmatrix(p)$Sigma, alpha_0 = 13, delta_0 = 8,
#' )
#' }
#' @export
gaussian_Gibbs <- function(X, y, b_0, B_0, alpha_0, delta_0,
                           beta, sigma, num_steps = 1e4, burn_ins = 1e3) {
  if (missing(beta))
    beta <- MASS::mvrnorm(1, b_0, B_0)
  if (missing(sigma))
    sigma <- 1 / sqrt(rgamma(1, alpha_0 / 2, delta_0 / 2))
  keep <- num_steps - burn_ins

  # Initialisation
  n <- length(y)
  alpha_1 <- alpha_0 + n
  beta_g <- beta
  sigma_g <- sigma
  inv_B_0 <- solve(B_0)
  inv_B_0_times_b_0 <- inv_B_0 %*% b_0
  XTX <- t(X) %*% X
  XTy <- t(X) %*% y

  res <- list()
  for (i in 1:num_steps) {
    # Update beta
    B_g <- solve(sigma_g^(-2) * XTX + inv_B_0)
    b_g <- B_g %*% (sigma_g^(-2) * XTy + inv_B_0_times_b_0)
    beta_g <- MASS::mvrnorm(1, b_g, B_g)

    # Update sigma
    delta_g <- delta_0 + sum((y - X %*% beta_g)^2)
    sigma_g <- 1 / sqrt(rgamma(1, alpha_1 / 2, delta_g / 2))

    # Keep track
    res <- push(res, list(beta = beta_g, sigma = sigma_g))
  }

  # Tidy format
  list(
    sigma = res %>% purrr::map_dbl(~.x$sigma) %>% tail(keep),
    beta = res %>% purrr::map(~.x$beta) %>% tail(keep) %>% do.call(rbind, .)
  )
}
