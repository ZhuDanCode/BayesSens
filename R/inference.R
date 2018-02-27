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
#' @param burn_ins Integer; number of burn-ins.
#' @examples
#' \dontrun{
#' n <- 100
#' p <- 10
#' data0 <- gaussian_data(n, p, intercept = T)
#' res <- gaussian_Gibbs(data0$X, data0$y,
#'   b_0 = rnorm(p+1), B_0 = pdmatrix(p+1)$Sigma,  # add one for the intercept
#'   alpha_0 = 13, delta_0 = 8,
#' )
#' }
#' @export
gaussian_Gibbs <- function(X, y, b_0, B_0, alpha_0, delta_0,
                           init_sigma, num_steps = 1e4, burn_ins = 1e3) {
  if (missing(init_sigma))
    init_sigma <- 1 / sqrt(rgamma(1, alpha_0 / 2, delta_0 / 2))

  # Initialisation
  n <- length(y)
  alpha_1 <- alpha_0 + n
  sigma_g <- init_sigma
  inv_B_0 <- solve(B_0)
  inv_B_0_times_b_0 <- inv_B_0 %*% b_0
  XTX <- t(X) %*% X
  XTy <- t(X) %*% y
  keep <- num_steps - burn_ins
  res <- vector("list", num_steps)

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
    res[[i]] <- list(beta = beta_g, sigma = sigma_g)
    setTxtProgressBar(pb, i)
  }

  # Tidy format
  list(
    sigma = res %>% purrr::map_dbl(~.x$sigma) %>% tail(keep),
    beta = res %>% purrr::map(~.x$beta) %>% tail(keep) %>% do.call(rbind, .)
  )
}


#' Model inference for student-t regression model with multivariate normal
#' prior for the mean and inverse gamma for the variance.
#' @param X A numeric matrix; the covariates.
#' @param y A numeric vector; the response variable.
#' @param b_0 A numeric vector; the mean for the multivariate normal prior.
#' @param B_0 A numeric matrix; the covariance for the multivariate normal prior.
#' @param alpha_0 A positive number; the shape parameter of the inverse-gamma prior.
#' @param delta_0 A positive number; the rate parameter of the inverse-gamma prior.
#' @param nu A positive number; the shape and rate parameters for the gamma prior.
#' @param init_beta (Optional) A numeric vector; the starting value of the beta parameter.
#' @param init_sigma (Optional) A positive number; the starting value of the sigma parameter.
#' @param num_steps Integer; number of MCMC steps.
#' @param burn_ins Integer; number of burn-ins.
#' @examples
#' \dontrun{
#' n <- 100
#' p <- 10
#' data0 <- student_t_data(n, p, intercept = TRUE)
#' res <- student_t_Gibbs(data0$X, data0$y,
#'   b_0 = rnorm(p+1), B_0 = pdmatrix(p+1)$Sigma,  # add one for intercept
#'   alpha_0 = 13, delta_0 = 8, nu = 5
#' )
#' }
#' @export
student_t_Gibbs <- function(X, y, b_0, B_0, alpha_0, delta_0, nu,
                    init_beta, init_sigma, num_steps = 1e4, burn_ins = 1e3) {
  if (missing(init_sigma))
    init_sigma <- 1 / sqrt(rgamma(1, alpha_0 / 2, delta_0 / 2))
  if (missing(init_beta))
    init_beta <- MASS::mvrnorm(1, b_0, B_0)

  # Initialisation
  n <- length(y)
  alpha_1 <- alpha_0 + n
  nu_1 <- nu + 1
  beta_g <- init_beta
  sigma_g <- init_sigma
  inv_B_0 <- solve(B_0)
  inv_B_0_times_b_0 <- inv_B_0 %*% b_0
  keep <- num_steps - burn_ins
  res <- vector("list", num_steps)

  pb <- txtProgressBar(1, num_steps, style = 3)
  for (i in 1:num_steps) {
    #Update lambda
    nu_2 <- nu + sigma_g^(-2) * (y - X %*% beta_g) ^ 2
    Lambda <- rgamma(n, nu_1 / 2, nu_2 / 2)

    # Update beta
    XTL <- left_multiply_D(t(X), Lambda)
    XTLX <- XTL %*% X
    XTLy <- XTL %*% y
    B_g <- solve(sigma_g^(-2) * XTLX + inv_B_0)
    b_g <- B_g %*% (sigma_g^(-2) * XTLy + inv_B_0_times_b_0)
    beta_g <- MASS::mvrnorm(1, b_g, B_g)

    # Update sigma
    delta_g <- delta_0 + left_multiply_D(t(y - X %*% beta_g), Lambda) %*% (y - X %*% beta_g)
    sigma_g <- 1 / sqrt(rgamma(1, alpha_1 / 2, delta_g / 2))

    # Keep track
    res[[i]] <- list(beta = beta_g, sigma = sigma_g)
    setTxtProgressBar(pb, i)
  }

  # Tidy format
  list(
    sigma = res %>% purrr::map_dbl(~.x$sigma) %>% tail(keep),
    beta = res %>% purrr::map(~.x$beta) %>% tail(keep) %>% do.call(rbind, .)
  )
}


