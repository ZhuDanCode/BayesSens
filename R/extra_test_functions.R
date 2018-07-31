# Sensitivity analysis for normal regression model with multivariate normal
# prior for the mean and inverse gamma for the variance.
#' @keywords internal
gaussian_AD_2 <- function(X, y, b_0, B_0, alpha_0, delta_0,
                        init_sigma, num_steps = 1e3) {
  if (missing(init_sigma))
    init_sigma <- 1 / sqrt(rgamma(1, alpha_0 / 2, delta_0 / 2))

  # Initialisation
  n <- length(y)
  len_beta <- length(b_0)
  alpha_1 <- alpha_0 + n
  sigma_g <- init_sigma
  runs_param <- vector("list", num_steps)

  # Autodiff helper variables and functions
  inv_B_0 <- solve(B_0)
  inv_B_0_times_b_0 <- inv_B_0 %*% b_0
  XTX <- crossprod(X)
  XTy <- crossprod(X, y)
  kp_B_0 <- t(inv_B_0) %x% inv_B_0
  d_transpose_K_nq <- commutation_matrix(length(y), 1)
  d_chol_I_n <- Matrix::Diagonal(nrow(B_0))
  d_chol_I_nn <- Matrix::Diagonal(nrow(B_0)^2)

  param_len <- c(length(b_0), length(B_0), length(alpha_0), length(delta_0),
                 length(init_sigma))
  init_differential <- function(len0) {
    param_len %>%
      purrr::map(~zeros(len0, .x)) %>%
      purrr::set_names(c("d_b0", "d_B0", "d_alpha0", "d_delta0", "d_sigma2_0"))
  }
  d_b0 <- init_differential(length(b_0))
  d_b0$d_b0 <- Matrix::Diagonal(length(b_0))
  d_B0 <- init_differential(length(B_0))
  d_B0$d_B0 <- Matrix::Diagonal(length(B_0))
  d_inv_B_0 <- d_inv(inv_B_0, d_B0)
  d_inv_B_0_times_b_0 <- d_product(inv_B_0, d_inv_B_0, b_0, d_b0)
  d_delta0 <- init_differential(length(delta_0))
  d_delta0$d_delta0 <- matrix(1)
  d_sigma2 <- init_differential(length(init_sigma))
  d_sigma2$d_sigma2_0 <- matrix(1)
  dG <- init_differential(1)
  d_XTX <- 0
  d_XTy <- 0
  d_X <- d_y <- d_z <- 0

  runs_d_beta <- vector("list", num_steps)
  runs_d_sigma2 <- vector("list", num_steps)

  # MCMC and Autodiff Loop
  pb <- txtProgressBar(1, num_steps, style = 3)
  for (i in 1:num_steps) {
    # Update beta
    sigma_g_n2 <- as.numeric(sigma_g^(-2))
    d_sigma_g_n2 <- d_inv(sigma_g_n2, d_sigma2)

    inv_sigma_XTX_p_inv_B_0 <- sigma_g_n2 * XTX + inv_B_0
    d_inv_sigma_XTX_p_inv_B_0 <- d_constant_multiply_matrix(
      sigma_g_n2, d_sigma_g_n2, XTX, d_XTX
    ) %>%
      d_sum(d_inv_B_0)

    B_g <- solve(inv_sigma_XTX_p_inv_B_0)
    d_Bg <- d_inv(B_g, d_inv_sigma_XTX_p_inv_B_0)

    b_g <- B_g %*% (sigma_g_n2 * XTy + inv_B_0_times_b_0)
    d_bg <- d_product(
      B_g, d_Bg,
      (sigma_g_n2 * XTy + inv_B_0_times_b_0),
      d_sum(
        d_constant_multiply_matrix(sigma_g_n2, d_sigma_g_n2, XTy, d_XTy),
        d_inv_B_0_times_b_0
      )
    )

    z <- rnorm(len_beta)
    chol_Bg <- t(chol(B_g))
    d_chol_Bg <- d_chol(chol_Bg, d_Bg, I_n = d_chol_I_n, I_nn = d_chol_I_nn)

    beta_g <- b_g + chol_Bg %*% z
    d_beta <- d_sum(d_bg, d_product(chol_Bg, d_chol_Bg, z, d_z))

    # Update sigma
    y_minus_X_beta_g <- y - X %*% beta_g
    d_y_minus_X_beta_g <- d_minus(d_y, d_product(X, d_X, beta_g, d_beta))
    d_t_y_minus_X_beta_g <- d_transpose(y_minus_X_beta_g, d_y_minus_X_beta_g, d_transpose_K_nq)

    delta_g <- as.numeric(delta_0 + crossprod(y_minus_X_beta_g))
    # delta_g <- delta_0 + t(y - X %*% beta_g) %*% (y - X %*% beta_g)
    d_delta_g <- d_sum(
      d_delta0,
      d_XXT(t(y_minus_X_beta_g), d_t_y_minus_X_beta_g)
    )

    G <- rgamma(1, alpha_1 / 2, 1)
    dG$d_alpha0 <- as.matrix(0.5 * d_Gamma(G, alpha_1 / 2))

    sigma_g <- 1 / sqrt(2 / delta_g * G)
    # sigma_g2 <- delta_g / (2 * G)
    d_sigma2 <- d_constant_multiply_constant(
      delta_g, d_delta_g,
      1 / (2*G), d_constant_inv(1 / (2*G), d_constant_multiply_constant(2, 0, G, dG))
    )

    # Keep track
    runs_param[[i]] <- list(beta = beta_g, sigma = sigma_g)
    runs_d_beta[[i]] <- d_beta
    runs_d_sigma2[[i]] <- d_sigma2
    setTxtProgressBar(pb, i)
  }

  # Tidy format
  append(
    tidy_list(runs_param),
    list(
      d_sigma2 = tidy_list(runs_d_sigma2),
      d_beta = tidy_list(runs_d_beta)
    )
  )
}


# Wishart Gibbs with Wishart draws implemented via Bartlett's decomposition.
#' @keywords internal
wishart_test_fun <- function(Xy, y, b_0, B_0, Xs, s, g_0, G_0, v_0, R_0,
                             init_gamma, init_sigma, num_steps = 1e3) {
  if (missing(init_gamma))
    init_gamma <- g_0 + t(chol(G_0)) %*% rnorm(length(g_0)) %>% as.vector()
  if (missing(init_sigma))
    init_sigma <- rWishart(1, v_0, R_0)

  # Initialisation
  n <- length(y)
  len_b0 <- length(b_0)
  len_g0 <- length(g_0)
  v_1 <- v_0 + n
  gamma_g <- init_gamma
  sigma_g <- init_sigma
  inv_B_0 <- solve(B_0)
  inv_B_0_times_b_0 <- inv_B_0 %*% b_0
  inv_G_0 <- solve(G_0)
  inv_G_0_times_g_0 <- inv_G_0 %*% g_0
  inv_R_0 <- solve(R_0)

  v <- cbind(y, s)
  X <- cbind(Xy, Xs)

  XyTXy <- crossprod(Xy)
  XyTXs <- crossprod(Xy, Xs)
  XsTXy <- crossprod(Xs, Xy)
  XsTXs <- crossprod(Xs)
  XyTy <- crossprod(Xy, y)
  XyTs <- crossprod(Xy, s)
  XsTy <- crossprod(Xs, y)
  XsTs <- crossprod(Xs, s)
  vTv <- crossprod(v, v)
  vTX <- crossprod(v, X)

  res <- vector("list", num_steps)

  #pre-simluation
  Zy <- matrix(rnorm(len_b0 * num_steps), len_b0, num_steps)
  Zs <- matrix(rnorm(len_g0 * num_steps), len_g0, num_steps)
  Z2 <- array(rnorm(2 * 2 * num_steps), dim = c(2, 2, num_steps)) #2 because there are 2 equations
  C <- matrix(0, 2, num_steps) #2 because there are 2 equations
  for (i in 1:2) {
    C[i,] <- rchisq(num_steps, v_1 - i + 1)
  }
  C <- sqrt(C)

  pb <- txtProgressBar(1, num_steps, style = 3)
  for (i in 1:num_steps) {
    #Update beta
    w_11_g <- sigma_g[1] - sigma_g[2] * sigma_g[3]/ sigma_g[4]
    b_sum <- XyTy - sigma_g[2] / sigma_g[4] * (XyTs - XyTXs %*% gamma_g)
    B_g <- solve(XyTXy / w_11_g + inv_B_0)
    b_g <- B_g %*% (b_sum / w_11_g + inv_B_0_times_b_0)
    chol_Bg <- t(chol(B_g))
    beta_g <- b_g + chol_Bg %*% Zy[, i]

    # Update gamma
    w_22_g <- sigma_g[4] - sigma_g[2] * sigma_g[3] / sigma_g[1]
    g_sum <- XsTs - sigma_g[2] / sigma_g[1] * (XsTy - XsTXy %*% beta_g)
    G_g <- solve(XsTXs / w_22_g + inv_G_0)
    g_g <- G_g %*% (g_sum / w_22_g + inv_G_0_times_g_0)
    chol_Gg <- t(chol(G_g))
    gamma_g <- g_g + chol_Gg %*% Zs[, i]

    # Update sigma
    delta <- matrix(c(beta_g, rep(0, len_b0 + len_g0), gamma_g), ncol = 2)
    vTX_delta <- vTX %*% delta
    X_delta <- X %*% delta
    R_1 <- solve(R_0 + vTv - vTX_delta - t(vTX_delta) + crossprod(X_delta))
    L <- t(chol(R_1))
    A <- diag(C[, i]) + get_lower_tri(Z2[,,i])
    LA <- L %*% A
    sigma_g <- solve(LA %*% t(LA))

    # Keep track
    res[[i]] <- list(beta = beta_g, gamma = gamma_g, sigma = sigma_g)
    setTxtProgressBar(pb, i)
  }

  # Tidy format
  list(
    sigma = res %>% purrr::map(~.x$sigma) %>% lapply(as.vector) %>% do.call(rbind, .),
    beta = res %>% purrr::map(~.x$beta) %>% lapply(as.vector) %>% do.call(rbind, .),
    gamma = res %>% purrr::map(~.x$gamma) %>% lapply(as.vector) %>% do.call(rbind, .)
  )
}


# Helper functions for testing
colMeans_tail <- function(X, p = 0.9) {
  colMeans(tail(X, max(round(nrow(X) * p), 1)))
}

mean_tail <- function(v0, p = 0.9) {
  mean(tail(v0, max(round(length(v0) * p), 1)))
}
