#' Sensitivity analysis for student-t regression model with multivariate normal
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
#' res <- student_t_AD(data0$X, data0$y,
#'   b_0 = rnorm(p+1), B_0 = pdmatrix(p+1)$Sigma,  # add one for intercept
#'   alpha_0 = 13, delta_0 = 8, nu = 5
#' )
#' }
#' @export
student_t_AD <- function(X, y, b_0, B_0, alpha_0, delta_0, nu,
                            init_beta, init_sigma, num_steps = 1e4, burn_ins = 1e3) {
  if (missing(init_sigma))
    init_sigma <- 1 / sqrt(rgamma(1, alpha_0 / 2, delta_0 / 2))
  if (missing(init_beta))
    init_beta <- MASS::mvrnorm(1, b_0, B_0)

  # Initialisation
  n <- length(y)
  len_beta <- length(init_beta)
  alpha_1 <- alpha_0 + n
  nu_1 <- nu + 1
  beta_g <- init_beta
  sigma_g <- init_sigma
  inv_B_0 <- solve(B_0)
  inv_B_0_times_b_0 <- inv_B_0 %*% b_0
  kp_B_0 <- t(inv_B_0) %x% inv_B_0
  keep <- num_steps - burn_ins
  assertthat::assert_that(
    keep > 0, msg = sprintf("num_steps = %d <= %d = burn_ins", num_steps, burn_ins)
  )

  d_beta <- init_student_t_differential(len_beta, len_beta)
  d_beta$d_beta0 <- diag(len_beta)
  d_sigma2 <- init_student_t_differential(1, len_beta)
  d_sigma2$d_sigma2_0 <- matrix(1)
  d_Lambda <- init_student_t_differential(n^2, len_beta)
  I_n <- diag(len_beta)
  I_nn <- diag(len_beta^2)
  K_nn <- matrixcalc::commutation.matrix(len_beta)
  elimL <- matrixcalc::elimination.matrix(len_beta)
  runs_param <- vector("list", num_steps)
  runs_d_beta <- vector("list", num_steps)
  runs_d_sigma2 <- vector("list", num_steps)

  # Auto-diff functions
  # Initialise storage variables externally to improve performance
  ext_d_Bg <- init_student_t_differential(len_beta^2, len_beta)
  ext_d_bg <- init_student_t_differential(len_beta, len_beta)
  ext_d_beta <- ext_d_bg
  deriv_Bg <- function(sigma_g, d_sigma2, Lambda_g, d_Lambda) {
    d_Bg <- ext_d_Bg
    XTLX <- A_times_diag_v0(t(X), as.numeric(Lambda_g)) %*% X
    inv_A <- solve(XTLX / sigma_g^2 + inv_B_0)  # intermediate variable
    fac_1 <- - (t(inv_A) %x% inv_A)
    fac_2 <- - matrixcalc::vec(XTLX) / sigma_g^4
    fac_2b <- sigma_g^(-2) * (t(X) %x% t(X))
    fac_3 <- fac_1 %*% fac_2
    fac_3b <- fac_1 %*% fac_2b

    d_Bg$d_b0 <- fac_3 %*% d_sigma2$d_b0 + fac_3b %*% d_Lambda$d_b0
    d_Bg$d_B0 <- fac_1 %*% (fac_2 %*% d_sigma2$d_B0 + fac_2b %*% d_Lambda$d_B0 - kp_B_0)
    d_Bg$d_alpha0 <- fac_3 %*% d_sigma2$d_alpha0 + fac_3b %*% d_Lambda$d_alpha0
    d_Bg$d_delta0 <- fac_3 %*% d_sigma2$d_delta0 + fac_3b %*% d_Lambda$d_delta0
    d_Bg$d_nu <- fac_3 %*% d_sigma2$d_nu + fac_3b %*% d_Lambda$d_nu
    d_Bg$d_beta0 <- fac_3 %*% d_sigma2$d_beta0 + fac_3b %*% d_Lambda$d_beta0
    d_Bg$d_sigma2_0 <- fac_3 %*% d_sigma2$d_sigma2_0 + fac_3b %*% d_Lambda$d_sigma2_0
    d_Bg
  }
  deriv_bg <- function(sigma_g, d_sigma2, B_g, d_Bg, Lambda_g, d_Lambda) {
    d_bg <- ext_d_bg
    XTLy <- A_times_diag_v0(t(X), as.numeric(Lambda_g)) %*% y
    fac_1 <- XTLy / sigma_g^2 + inv_B_0_times_b_0
    tfac_1 <- (t(fac_1) %x% diag(nrow(B_g)))
    fac_2 <- - (XTLy) / sigma_g^4
    fac_3 <- t(y) %x% t(X) / sigma_g^2
    d_bg$d_b0 <- tfac_1 %*% d_Bg$d_b0 + B_g %*% (fac_2 %*% d_sigma2$d_b0 + fac_3 %*% d_Lambda$d_b0 + inv_B_0)
    d_bg$d_B0 <- tfac_1 %*% d_Bg$d_B0 + B_g %*%
      (fac_2 %*% d_sigma2$d_B0 + fac_3 %*% d_Lambda$d_B0 - t(inv_B_0_times_b_0) %x% inv_B_0)
    d_bg$d_alpha0 <- tfac_1 %*% d_Bg$d_alpha0 + B_g %*% (fac_2 %*% d_sigma2$d_alpha0 + fac_3 %*% d_Lambda$d_alpha0)
    d_bg$d_delta0 <- tfac_1 %*% d_Bg$d_delta0 + B_g %*% (fac_2 %*% d_sigma2$d_delta0 + fac_3 %*% d_Lambda$d_delta0)
    d_bg$d_nu <- tfac_1 %*% d_Bg$d_nu + B_g %*% (fac_2 %*% d_sigma2$d_nu + fac_3 %*% d_Lambda$d_nu)
    d_bg$d_beta0 <- tfac_1 %*% d_Bg$d_beta0 + B_g %*% (fac_2 %*% d_sigma2$d_beta0 + fac_3 %*% d_Lambda$d_beta0)
    d_bg$d_sigma2_0 <- tfac_1 %*% d_Bg$d_sigma2_0 + B_g %*% (fac_2 %*% d_sigma2$d_sigma2_0  + fac_3 %*% d_Lambda$d_sigma2_0)
    d_bg
  }
  deriv_beta <- function(sigma_g, d_sigma2, Lambda_g, d_Lambda, B_g, z) {
    n <- nrow(B_g)
    L <- t(chol(B_g))
    fac_1 <- t(elimL) %*% solve(elimL %*% (I_nn + K_nn) %*% (L %x% I_n) %*% t(elimL)) %*% elimL
    # fac_1 <- solve((L %x% I_n) + (I_n %x% L) %*% K_nn)
    z_mat <- t(z) %x% I_n %*% fac_1

    d_beta <- ext_d_beta
    d_Bg <- deriv_Bg(sigma_g, d_sigma2, Lambda_g, d_Lambda)
    d_bg <- deriv_bg(sigma_g, d_sigma2, B_g, d_Bg, Lambda_g, d_Lambda)

    d_beta$d_b0 <- d_bg$d_b0 +  z_mat %*% d_Bg$d_b0
    d_beta$d_B0 <- d_bg$d_B0 + z_mat %*% d_Bg$d_B0
    d_beta$d_alpha0 <- d_bg$d_alpha0 + z_mat %*% d_Bg$d_alpha0
    d_beta$d_delta0 <- d_bg$d_delta0 + z_mat %*% d_Bg$d_delta0
    d_beta$d_nu <- d_bg$d_nu + z_mat %*% d_Bg$d_nu
    d_beta$d_beta0 <- d_bg$d_beta0 + z_mat %*% d_Bg$d_beta0
    d_beta$d_sigma2_0 <- d_bg$d_sigma2_0 + z_mat %*% d_Bg$d_sigma2_0
    d_beta
  }

  ext_d_delta <- init_student_t_differential(1, len_beta)
  ext_d_sigma2 <- ext_d_G <- ext_d_delta
  deriv_delta <- function(beta_g, d_beta, Lambda_g, d_Lambda) {
    d_delta <- ext_d_delta
    fac_1 <- - 2 * A_times_diag_v0(t(y - X %*% beta_g), Lambda_g) %*% X
    fac_2 <- (y - X %*% beta_g)
    fac_2 <- t(fac_2) %x% t(fac_2)

    d_delta$d_b0 <- fac_1 %*% d_beta$d_b0 + fac_2 %*% d_Lambda$d_b0
    d_delta$d_B0 <- fac_1 %*% d_beta$d_B0 + fac_2 %*% d_Lambda$d_B0
    d_delta$d_alpha0 <- fac_1 %*% d_beta$d_alpha0 + fac_2 %*% d_Lambda$d_alpha0
    d_delta$d_delta0 <- 1 + fac_1 %*% d_beta$d_delta0 + fac_2 %*% d_Lambda$d_delta0
    d_delta$d_nu <- fac_1 %*% d_beta$d_nu + fac_2 %*% d_Lambda$d_nu
    d_delta$d_beta0 <- fac_1 %*% d_beta$d_beta0 + fac_2 %*% d_Lambda$d_beta0
    d_delta$d_sigma2_0 <- fac_1 %*% d_beta$d_sigma2_0 + fac_2 %*% d_Lambda$d_sigma2_0
    d_delta
  }
  deriv_G <- function(G, alpha, var_name) {
    d_G <- ext_d_G
    f <- function(t) { log(t) * dgamma(t, alpha, 1) }
    num_1 <- integrate(f, 0, G)$value
    num_2 <- digamma(alpha) * pgamma(G, alpha, 1)
    d_G[[var_name]] <- - 0.5 * (num_1 - num_2) / dgamma(G, alpha, 1)
    d_G
  }
  deriv_sigma2 <- function(beta_g, delta_g, G, alpha_1, d_beta, Lambda_g, d_Lambda) {
    delta_g <- as.numeric(delta_g)
    d_delta <- deriv_delta(beta_g, d_beta, Lambda_g, d_Lambda)
    d_G <- deriv_G(G, 0.5 * alpha_1, "d_alpha0")
    d_sigma2 <- ext_d_sigma2
    d_sigma2$d_b0 <- d_delta$d_b0 / (2 * G) - delta_g / (2 * G^2) * d_G$d_b0
    d_sigma2$d_B0 <- d_delta$d_B0 / (2 * G) - delta_g / (2 * G^2) * d_G$d_B0
    d_sigma2$d_alpha0 <- d_delta$d_alpha0 / (2 * G) - delta_g / (2 * G^2) * d_G$d_alpha0
    d_sigma2$d_delta0 <- d_delta$d_delta0 / (2 * G) - delta_g / (2 * G^2) * d_G$d_delta0
    d_sigma2$d_nu <- d_delta$d_nu / (2 * G) - delta_g / (2 * G^2) * d_G$d_nu
    d_sigma2$d_beta0 <- d_delta$d_beta0 / (2 * G) - delta_g / (2 * G^2) * d_G$d_beta0
    d_sigma2$d_sigma2_0 <- d_delta$d_sigma2_0 / (2 * G) - delta_g / (2 * G^2) * d_G$d_sigma2_0
    d_sigma2
  }

  ext_d_v_2i <- ext_d_lambda_i <- ext_d_delta
  ext_d_Lambda <- init_student_t_differential(n^2, len_beta)
  deriv_v_2i <- function(i, nu, sigma_g, d_sigma2, beta_g, d_beta) {
    d_v_2i <- ext_d_v_2i
    fac_2 <- as.numeric(sigma_g^(-2) * (y[i] - X[i,] %*% beta_g))
    fac_1 <- as.numeric(fac_2^2)
    d_v_2i$d_b0 <- - fac_1 * d_sigma2$d_b0 - 2 * fac_2 * X[i,] %*% d_beta$d_b0
    d_v_2i$d_B0 <- - fac_1 * d_sigma2$d_B0 - 2 * fac_2 * X[i,] %*% d_beta$d_B0
    d_v_2i$d_alpha0 <- - fac_1 * d_sigma2$d_alpha0 - 2 * fac_2 * X[i,] %*% d_beta$d_alpha0
    d_v_2i$d_delta0 <- - fac_1 * d_sigma2$d_delta0 - 2 * fac_2 * X[i,] %*% d_beta$d_delta0
    d_v_2i$d_nu <- 1 - fac_1 * d_sigma2$d_nu - 2 * fac_2 * X[i,] %*% d_beta$d_nu
    d_v_2i$d_beta0 <- - fac_1 * d_sigma2$d_beta0 - 2 * fac_2 * X[i,] %*% d_beta$d_beta0
    d_v_2i$d_sigma2_0 <- - fac_1 * d_sigma2$d_sigma2_0 - 2 * fac_2 * X[i,] %*% d_beta$d_sigma2_0
    d_v_2i
  }
  deriv_lambda_i <- function(i, G, dG, nu, d_v_2i, sigma_g, d_sigma2, beta_g, d_beta) {
    d_lambda_i <- ext_d_lambda_i
    v_2i <- nu + sigma_g^(-2) * (y[i] - X[i,] %*% beta_g)^2
    fac1 <- as.numeric(- 2 / v_2i^2 * G)
    fac2 <- as.numeric(2 / v_2i)
    d_lambda_i$d_b0 <- fac1 * d_v_2i$d_b0 + fac2 * dG$d_b0
    d_lambda_i$d_B0 <- fac1 * d_v_2i$d_B0 + fac2 * dG$d_B0
    d_lambda_i$d_alpha0 <- fac1 * d_v_2i$d_alpha0 + fac2 * dG$d_alpha0
    d_lambda_i$d_delta0 <- fac1 * d_v_2i$d_delta0 + fac2 * dG$d_delta0
    d_lambda_i$d_nu <- fac1 * d_v_2i$d_nu + fac2 * dG$d_nu
    d_lambda_i$d_beta0 <- fac1 * d_v_2i$d_beta0 + fac2 * dG$d_beta0
    d_lambda_i$d_sigma2_0 <- fac1 * d_v_2i$d_sigma2_0 + fac2 * dG$d_sigma2_0
    d_lambda_i
  }
  deriv_Lambda <- function(nu, G, sigma_g, d_sigma2, beta_g, d_beta) {
    d_v_2i <- purrr::map(seq(n), ~deriv_v_2i(.x, nu, sigma_g, d_sigma2, beta_g, d_beta))
    dG <- purrr::map(seq(n), ~ deriv_G(G[.x], 0.5 * nu_1, "d_nu"))
    d_lambda_i <- seq(n) %>%
      purrr::map(~deriv_lambda_i(.x, G[.x], dG[[.x]], nu, d_v_2i[[.x]], sigma_g, d_sigma2, beta_g, d_beta)) %>%
      collect()
    d_Lambda <- ext_d_Lambda
    d_Lambda$d_b0 <- d_lambda_i$d_b0 %>% add_zeros()
    d_Lambda$d_B0 <- d_lambda_i$d_B0 %>% add_zeros()
    d_Lambda$d_alpha0 <- d_lambda_i$d_alpha0 %>% add_zeros()
    d_Lambda$d_delta0 <- d_lambda_i$d_delta0 %>% add_zeros()
    d_Lambda$d_nu <- d_lambda_i$d_nu %>% add_zeros()
    d_Lambda$d_beta0 <- d_lambda_i$d_beta0 %>% add_zeros()
    d_Lambda$d_sigma2_0 <- d_lambda_i$d_sigma2_0 %>% add_zeros()
    d_Lambda
  }
  add_zeros <- function(m0) {
    total_row <- (nrow(m0) - 1) * n + nrow(m0)
    m1 <- matrix(0, nrow = total_row, ncol = ncol(m0))
    m1[seq(1, total_row, n+1), ] <- m0
    m1
  }

  pb <- txtProgressBar(1, num_steps, style = 3)
  for (i in 1:num_steps) {
    #Update lambda
    nu_2 <- nu + sigma_g^(-2) * (y - X %*% beta_g) ^ 2
    g2 <- rgamma(n, nu_1 / 2, 1)
    Lambda_g <- 2 / nu_2 * g2
    #Autodiff gamma
    d_Lambda <- deriv_Lambda(nu, g2, sigma_g, d_sigma2, beta_g, d_beta)

    # Update beta
    XTL <- A_times_diag_v0(t(X), as.numeric(Lambda_g))
    XTLX <- XTL %*% X
    XTLy <- XTL %*% y
    B_g <- solve(sigma_g^(-2) * XTLX + inv_B_0)
    b_g <- B_g %*% (sigma_g^(-2) * XTLy + inv_B_0_times_b_0)
    z <- rnorm(len_beta)
    beta_g <- b_g + t(chol(B_g)) %*% z # originally, beta_g <- MASS::mvrnorm(1, b_g, B_g)
    # Autodiff beta
    d_beta <- deriv_beta(sigma_g, d_sigma2, Lambda_g, d_Lambda, B_g, z)

    # Update sigma
    delta_g <- as.numeric(delta_0 + A_times_diag_v0(t(y - X %*% beta_g), Lambda_g) %*% (y - X %*% beta_g))
    g <- rgamma(1, alpha_1 / 2, 1)
    sigma_g <- 1 / sqrt(2 / delta_g * g)
    # Autodiff sigma
    d_sigma2 <- deriv_sigma2(beta_g, delta_g, g, alpha_1, d_beta, Lambda_g, d_Lambda)

    # Keep track
    runs_param[[i]] <- list(beta = beta_g, sigma = sigma_g)
    runs_d_beta[[i]] <- d_beta
    runs_d_sigma2[[i]] <- d_sigma2
    setTxtProgressBar(pb, i)
  }

  # Tidy format
  list(
    sigma = runs_param %>% purrr::map_dbl(~.x$sigma) %>% tail(keep),
    beta = runs_param %>% purrr::map(~t(.x$beta)) %>% tail(keep) %>% do.call(rbind, .),
    d_sigma2 = runs_d_sigma2 %>% tail(keep) %>% tidy_student_t_differential(),
    d_beta = runs_d_beta %>% tail(keep) %>% tidy_student_t_differential()
  )
}


init_student_t_differential <- function(len0, l_b0) {
  list(
    d_b0 = zeros(len0, l_b0),
    d_B0 = zeros(len0, l_b0^2),
    d_alpha0 = zeros(len0, 1),
    d_delta0 = zeros(len0, 1),
    d_nu = zeros(len0, 1),
    d_beta0 = zeros(len0, l_b0),
    d_sigma2_0 = zeros(len0, 1)
  )
}


tidy_student_t_differential <- function(dlist0) {
  extract_rbind <- function(attr0) {
    dlist0 %>%
      purrr::map(~t(matrixcalc::vec(.x[[attr0]]))) %>%
      do.call(rbind, .)
  }
  list(
    d_b0 = extract_rbind("d_b0"),
    d_B0 = extract_rbind("d_B0"),
    d_alpha0 = extract_rbind("d_alpha0"),
    d_delta0 = extract_rbind("d_delta0"),
    d_nu = extract_rbind("d_nu"),
    d_beta0 = extract_rbind("d_beta0"),
    d_sigma2_0 = extract_rbind("d_sigma2_0")
  )
}
