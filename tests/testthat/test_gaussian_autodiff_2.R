context("Cross-checking the gaussian autodiff using a different implementation")
library(BayesSense)

gaussian_AD_2 <- function(X, y, b_0, B_0, alpha_0, delta_0, init_sigma, num_steps = 1e3) {
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

testthat::test_that("Test gaussian autodiff", {
  # skip_on_cran()
  skip_if_not(Sys.getenv("run_all_tests", TRUE))
  set.seed(123)
  n <- 1000
  p <- 5
  data0 <- gaussian_data(n, p, intercept = TRUE)

  set.seed(123)
  default_parameters <- list(
    X = data0$X, y = data0$y,
    b_0 = rnorm(p+1), B_0 = pdmatrix(p+1)$Sigma,  # add one for the intercept
    alpha_0 = 13, delta_0 = 8, num_steps = 200
  )

  set.seed(123)
  res <- do.call(gaussian_AD, default_parameters)
  set.seed(123)
  res2 <- do.call(gaussian_AD_2, default_parameters)

  for (i in names(res$d_beta)) {
    v1 <- colMeans_tail(res$d_beta[[i]])
    v2 <- colMeans_tail(res2$d_beta[[i]])
    expect_lt(sum(abs(v1 - v2)), 1e-10)
    print(rbind(implementation_1 = v1, implementation_2 = v2))
  }

  for (i in names(res$d_sigma2)) {
    v1 <- colMeans_tail(res$d_sigma2[[i]])
    v2 <- colMeans_tail(res2$d_sigma2[[i]])
    expect_lt(sum(abs(v1 - v2)), 1e-10)
    print(rbind(implementation_1 = v1, implementation_2 = v2))
  }
})
