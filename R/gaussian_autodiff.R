#' Sensitivity analysis for normal regression model with multivariate normal
#' prior for the mean and inverse gamma for the variance.
#' @param X A numeric matrix; the covariates.
#' @param y A numeric vector; the response variable.
#' @param b_0 A numeric vector; the mean for the multivariate normal prior.
#' @param B_0 A numeric matrix; the covariance for the multivariate normal prior.
#' @param alpha_0 A positive number; the shape parameter of the inverse-gamma prior.
#' @param delta_0 A positive number; the rate parameter of the inverse-gamma prior.
#' @param init_sigma (Optional) A numeric matrix
#' @param num_steps Integer; number of MCMC steps.
#' @examples
#' \dontrun{
#' n <- 1000
#' p <- 5
#' data0 <- gaussian_data(n, p, intercept = TRUE)
#' res <- gaussian_AD(
#'   X = data0$X, y = data0$y,
#'   b_0 = rnorm(p+1), B_0 = diag(p+1),  # add one for the intercept
#'   alpha_0 = 13, delta_0 = 8,
#' )
#' }
#' @export
gaussian_AD <- function(X, y, b_0, B_0, alpha_0, delta_0,
                        init_sigma, num_steps = 1e3) {
  if (missing(init_sigma))
    init_sigma <- 1 / sqrt(rgamma(1, alpha_0 / 2, delta_0 / 2))

  # Initialisation
  init_gauss_differential <- purrr::partial(
    init_differential,
    vec0 = purrr::map_dbl(list(b_0, B_0, alpha_0, delta_0, init_sigma), length),
    names0 = c("b0", "B0", "alpha0", "delta0", "sigma2_0")
  )
  n <- length(y)
  alpha_1 <- alpha_0 + n
  sigma_g <- init_sigma
  # Autodiff variables
  len_beta <- length(b_0)
  d_beta <- init_gauss_differential(length(b_0))
  d_sigma2 <- init_gauss_differential(length(sigma_g))
  d_sigma2$d_sigma2_0 <- matrix(1)
  # Variables to keep track
  runs_beta <- vector("list", num_steps)
  runs_sigma <- vector("list", num_steps)
  runs_d_beta <- vector("list", num_steps)
  runs_d_sigma2 <- vector("list", num_steps)

  # Helper variables and functions
  inv_B_0 <- solve(B_0)
  inv_B_0_times_b_0 <- inv_B_0 %*% b_0
  XTX <- crossprod(X)
  XTy <- crossprod(X, y)
  kp_B_0 <- t(inv_B_0) %x% inv_B_0
  I_n <- diag(len_beta)
  I_nn <- diag(len_beta^2)
  K_nn <- commutation_matrix(len_beta, len_beta)
  elimL <- elimination_matrix(len_beta)
  t_matrix <- . %>% as.matrix() %>% t()

  deriv_Bg <- function(sigma_g, d_sigma2) {
    d_Bg <- init_gauss_differential(length(B_0))
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
    d_bg <- init_gauss_differential(length(b_0))
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
    fac_1 <- t_matrix(elimL) %*% solve(elimL %*% (I_nn + K_nn) %*% (L %x% I_n) %*% t_matrix(elimL)) %*% elimL
    z_mat <- (t(z) %x% I_n) %*% fac_1

    d_beta <- init_gauss_differential(length(b_0))
    d_Bg <- deriv_Bg(sigma_g, d_sigma2)
    d_bg <- deriv_bg(sigma_g, d_sigma2, B_g, d_Bg)

    d_beta$d_b0 <- d_bg$d_b0 +  z_mat %*% d_Bg$d_b0
    d_beta$d_B0 <- d_bg$d_B0 + z_mat %*% d_Bg$d_B0
    d_beta$d_alpha0 <- d_bg$d_alpha0 + z_mat %*% d_Bg$d_alpha0
    d_beta$d_delta0 <- d_bg$d_delta0 + z_mat %*% d_Bg$d_delta0
    d_beta$d_sigma2_0 <- d_bg$d_sigma2_0 + z_mat %*% d_Bg$d_sigma2_0
    d_beta
  }
  deriv_delta <- function(beta_g, d_beta) {
    d_delta <- init_gauss_differential(length(delta_0))
    fac_1 <- - 2 * t(y - X %*% beta_g) %*% X
    d_delta$d_b0 <- fac_1 %*% d_beta$d_b0
    d_delta$d_B0 <- fac_1 %*% d_beta$d_B0
    d_delta$d_alpha0 <- fac_1 %*% d_beta$d_alpha0
    d_delta$d_delta0 <- 1 + fac_1 %*% d_beta$d_delta0
    d_delta$d_sigma2_0 <- fac_1 %*% d_beta$d_sigma2_0
    d_delta
  }
  deriv_G <- function(G, alpha) {
    d_G <- init_gauss_differential(1)
    f <- function(t) { log(t) * dgamma(t, alpha, 1) }
    num_1 <- integrate(f, 0, G)$value
    num_2 <- digamma(alpha) * pgamma(G, alpha, 1)
    d_G$d_alpha0 <- - 0.5 * (num_1 - num_2) / dgamma(G, alpha, 1)
    d_G
  }
  deriv_sigma2 <- function(beta_g, delta_g, G, alpha_1, d_beta) {
    d_delta <- deriv_delta(beta_g, d_beta)
    d_G <- deriv_G(G, 0.5 * alpha_1)
    d_sigma2 <- init_gauss_differential(length(sigma_g))
    d_sigma2$d_b0 <- d_delta$d_b0 / (2 * G) - delta_g / (2 * G^2) * d_G$d_b0
    d_sigma2$d_B0 <- d_delta$d_B0 / (2 * G) - delta_g / (2 * G^2) * d_G$d_B0
    d_sigma2$d_alpha0 <- d_delta$d_alpha0 / (2 * G) - delta_g / (2 * G^2) * d_G$d_alpha0
    d_sigma2$d_delta0 <- d_delta$d_delta0 / (2 * G) - delta_g / (2 * G^2) * d_G$d_delta0
    d_sigma2$d_sigma2_0 <- d_delta$d_sigma2_0 / (2 * G) - delta_g / (2 * G^2) * d_G$d_sigma2_0
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
    runs_beta[[i]] <- beta_g
    runs_sigma[[i]] <- sigma_g
    runs_d_beta[[i]] <- d_beta
    runs_d_sigma2[[i]] <- d_sigma2
    setTxtProgressBar(pb, i)
  }

  # Tidy format
  list(
    beta = map_reduce(runs_beta, t, rbind),
    sigma = map_reduce(runs_sigma, t, rbind),
    d_beta = tidy_list(runs_d_beta),
    d_sigma2 = tidy_list(runs_d_sigma2)
  )
}
