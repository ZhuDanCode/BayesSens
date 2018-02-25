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
#' data0 <- gaussian_data(n, p, intercept = T)
#' res <- gaussian_AD(data0$X, data0$y,
#'   b_0 = rnorm(p+1), B_0 = pdmatrix(p+1)$Sigma,  # add one for the intercept
#'   alpha_0 = 13, delta_0 = 8,
#' )
#' }
#' @export
gaussian_AD <- function(X, y, b_0, B_0, alpha_0, delta_0,
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
  # Autodiff variables
  len_beta <- length(beta)
  d_beta <- init_differential(len_beta, len_beta)
  d_sigma2 <- init_differential(1, len_beta)

  # Variables to keep track
  runs_param <- list()
  runs_d_beta <- list()
  runs_d_sigma2 <- list()
  # Helper variables and functions
  inv_B_0 <- solve(B_0)
  inv_B_0_times_b_0 <- inv_B_0 %*% b_0
  XTX <- t(X) %*% X
  XTy <- t(X) %*% y
  kp_B_0 <- t(inv_B_0) %x% inv_B_0
  deriv_Bg <- function(sigma_g, d_sigma2) {
    d_Bg <- init_differential(len_beta^2, len_beta)
    inv_A <- solve(XTX / sigma_g^2 + inv_B_0)  # intermediate variable
    fac_1 <- - (t(inv_A) %x% inv_A)
    fac_2 <- - matrixcalc::vec(XTX) / sigma_g^4
    fac_3 <- fac_1 %*% fac_2

    d_Bg$d_b0 <- fac_3 %*% d_sigma2$d_b0
    d_Bg$d_B0 <- fac_1 %*% (fac_2 %*% d_sigma2$d_B0 - kp_B_0)
    d_Bg$d_alpha0 <- fac_3 %*% d_sigma2$d_alpha0
    d_Bg$d_delta0 <- fac_3 %*% d_sigma2$d_delta0
    d_Bg$d_sigma2_0 <- fac_3 %*% d_sigma2$d_sigma2_0
    d_Bg
  }
  deriv_bg <- function(sigma_g, d_sigma2, B_g, d_Bg) {
    d_bg <- init_differential(len_beta, len_beta)
    fac_1 <- XTy / sigma_g^2 + inv_B_0_times_b_0
    tfac_1 <- (t(fac_1) %x% diag(nrow(B_g)))
    fac_2 <- - (XTy) / sigma_g^4
    d_bg$d_b0 <- tfac_1 %*% d_Bg$d_b0 + B_g %*% (fac_2 %*% d_sigma2$d_b0 + inv_B_0)
    d_bg$d_B0 <- tfac_1 %*% d_Bg$d_B0 + B_g %*%
      (fac_2 %*% d_sigma2$d_B0 - t(inv_B_0_times_b_0) %x% inv_B_0)
    d_bg$d_alpha0 <- tfac_1 %*% d_Bg$d_alpha0 + B_g %*% (fac_2 %*% d_sigma2$d_alpha0)
    d_bg$d_delta0 <- tfac_1 %*% d_Bg$d_delta0 + B_g %*% (fac_2 %*% d_sigma2$d_delta0)
    d_bg$d_sigma2_0 <- tfac_1 %*% d_Bg$d_sigma2_0 + B_g %*% (fac_2 %*% d_sigma2$d_sigma2_0)
    d_bg
  }
  deriv_beta <- function(sigma_g, d_sigma2, B_g, z) {
    n <- nrow(B_g)
    L <- t(chol(B_g))
    I_n <- diag(n)
    I_nn <- diag(n^2)
    K_nn <- matrixcalc::commutation.matrix(n)
    elimL <- matrixcalc::elimination.matrix(n)
    commD <- matrixcalc::duplication.matrix(n)
    fac_1 <- t(elimL) %*% solve(elimL %*% (I_nn + K_nn) %*% (L %x% I_n) %*% commD) %*% elimL
    # fac_1 <- solve((L %x% I_n) + (I_n %x% L) %*% K_nn)
    z_mat <- t(z) %x% I_n

    d_beta <- init_differential(len_beta, len_beta)
    d_Bg <- deriv_Bg(sigma_g, d_sigma2)
    d_bg <- deriv_bg(sigma_g, d_sigma2, B_g, d_Bg)

    d_beta$d_b0 <- d_bg$d_b0 +  z_mat %*% fac_1 %*% d_Bg$d_b0
    d_beta$d_B0 <- d_bg$d_B0 + z_mat %*% fac_1 %*% d_Bg$d_B0
    d_beta$d_alpha0 <- d_bg$d_alpha0 + z_mat %*% fac_1 %*% d_Bg$d_alpha0
    d_beta$d_delta0 <- d_bg$d_delta0 + z_mat %*% fac_1 %*% d_Bg$d_delta0
    d_beta$d_sigma2_0 <- d_bg$d_sigma2_0 + z_mat %*% fac_1 %*% d_Bg$d_sigma2_0
    d_beta
  }
  deriv_delta <- function(beta_g, d_beta) {
    d_delta <- init_differential(1, len_beta)
    fac_1 <- - 2 * t(y - X %*% beta_g) %*% X
    d_delta$d_b0 <- fac_1 %*% d_beta$d_b0
    d_delta$d_B0 <- fac_1 %*% d_beta$d_B0
    d_delta$d_alpha0 <- fac_1 %*% d_beta$d_alpha0
    d_delta$d_delta0 <- 1 + fac_1 %*% d_beta$d_delta0
    d_delta$d_sigma2_0 <- fac_1 %*% d_beta$d_sigma2_0
    d_delta
  }
  deriv_G <- function(G, alpha) {
    d_G <- init_differential(1, len_beta)
    f <- function(t) { log(t) * t^(alpha - 1) * exp(-t) }
    integral_1 <- integrate(f, 0, G)
    num_1 <- integral_1$value / gamma(alpha)  #divide gamma(alpha) after to avoid repeated evaluation
    num_2 <- digamma(alpha) * pgamma(G, alpha, 1)
    d_G$d_alpha0 <- - (num_1 - num_2) / dgamma(G, alpha, 1)
    d_G
  }
  deriv_sigma2 <- function(beta_g, delta_g, G, alpha_1, d_beta) {
    d_delta <- deriv_delta(beta_g, d_beta)
    d_G <- deriv_G(G, 0.5 * alpha_1)
    d_sigma2 <- init_differential(1, len_beta)
    d_sigma2$d_b0 <- d_delta$d_b0 / (2 * G) - delta_g / (2 * G) * d_G$d_b0
    d_sigma2$d_B0 <- d_delta$d_B0 / (2 * G) - delta_g / (2 * G) * d_G$d_B0
    d_sigma2$d_alpha0 <- d_delta$d_alpha0 / (2 * G) - delta_g / (2 * G) * d_G$d_alpha0
    d_sigma2$d_delta0 <- d_delta$d_delta0 / (2 * G) - delta_g / (2 * G) * d_G$d_delta0
    d_sigma2$d_sigma2_0 <- d_delta$d_sigma2_0 / (2 * G) - delta_g / (2 * G) * d_G$d_sigma2_0
    d_sigma2
  }

  # MCMC and Autodiff Loop
  pb <- txtProgressBar(1, num_steps, style = 3)
  for (i in 1:num_steps) {
    # Update beta
    B_g <- solve(sigma_g^(-2) * XTX + inv_B_0)
    b_g <- B_g %*% (sigma_g^(-2) * XTy + inv_B_0_times_b_0)
    z <- rnorm(len_beta)
    beta_g <- b_g + t(chol(B_g)) %*% z # originally, beta_g <- MASS::mvrnorm(1, b_g, B_g)
    # Autodiff beta
    d_beta <- deriv_beta(sigma_g, d_sigma2, B_g, z)

    # Update sigma
    delta_g <- delta_0 + sum((y - X %*% beta_g)^2)
    G <- rgamma(1, alpha_1 / 2, 1)
    sigma_g <- 1 / sqrt(2 / delta_g * G)
    # Autodiff sigma
    d_sigma2 <- deriv_sigma2(beta_g, delta_g, G, alpha_1, d_beta)

    # Keep track
    runs_param %<>% push(list(beta = beta_g, sigma = sigma_g))
    runs_d_beta %<>% push(d_beta)
    runs_d_sigma2 %<>% push(d_sigma2)
    setTxtProgressBar(pb, i)
  }

  # Tidy format
  list(
    sigma = runs_param %>% purrr::map_dbl(~.x$sigma) %>% tail(keep),
    beta = runs_param %>% purrr::map(~.x$beta) %>% tail(keep) %>% do.call(rbind, .),
    d_sigma2 = runs_d_sigma2,
    d_beta = runs_d_beta
  )
}
