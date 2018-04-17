#' Sensitivity analysis for Vector-Auto-Regressive model with independent Normal-Wishart priors.
#' @param data0 Matrix; each row is one observation, each column is one measurement / predictor.
#' @param lag Integer; the lag of the time series model.
#' @param b_0 A numeric vector; the mean for the multivariate normal prior.
#' @param B_0 A numeric matrix; the covariance for the multivariate normal prior.
#' @param v_0 scalar; degree of freedom for the inverse-Wishart prior.
#' @param S_0 matrix; the scale matrix for the inverse-Wishart prior.
#' @param init_Sigma (Optional) matrix; the starting value of the noise inverse-covariance.
#' @param num_steps integer; number of MCMC steps.
#' @export
VAR_AD <- function(data0, lag, b_0, B_0, v_0, S_0, init_Sigma, num_steps = 3e3) {
  data0 <- train_data(data0, lag)
  Y <- data0$Y
  X <- data0$X

  T0 <- nrow(Y)
  n <- ncol(Y)
  p <- lag

  K <- 1 + n * p
  d <- n * K

  invB <- solve(B_0)
  invVb <- invB %*% b_0
  v_1 <- v_0 + T0
  if (missing(init_Sigma))
    init_Sigma <- MCMCpack::rwish(v_0, S_0)
  invSigma <- init_Sigma

  XTX <- t(X) %*% X
  XTY <- t(X) %*% Y
  YTY <- t(Y) %*% Y

  total <- num_steps
  Z <- matrix(rnorm(d * total), d, total)
  Z2 <- array(rnorm(n * n * total), dim = c(n, n, total))
  C <- matrix(0, n, total)
  for (i in 1:n) {
    C[i,] <- rchisq(total, v_1 - i + 1)
  }
  C <- sqrt(C)
  res_beta <- matrix(0, total, length(b_0))
  res_sigma <- matrix(0, total, length(init_Sigma))

  get_lower_tri <- . %>% {.[upper.tri(., T)] <- 0; .}

  # Autodiff
  init_deriv <- . %>% {list(
    d_b0 = zeros(., length(b_0)),
    d_B0 = zeros(., length(B_0)),
    d_nu0 = zeros(., length(v_0)),
    d_S0 = zeros(., length(S_0)),
    d_Sigma0 = zeros(., length(init_Sigma))
  )}
  d_b0 <- init_deriv(length(b_0))
  diag(d_b0$d_b0) <- 1
  d_B0 <- init_deriv(length(B_0))
  diag(d_B0$d_B0) <- 1
  d_S0 <- init_deriv(length(S_0))
  diag(d_S0$d_S0) <- 1
  d_invSigma <- init_deriv(length(init_Sigma))
  d_invSigma$d_Sigma0 <- t(invSigma) %x% invSigma
  d_XTX <- 0
  d_XTY <- 0
  d_YTY <- 0
  d_Z_i <- 0
  runs_d_beta <- vector("list", total)
  runs_d_Sigma <- vector("list", total)
  I_m <- Matrix::Diagonal(nrow(B_0))
  I_mm <- Matrix::Diagonal(nrow(B_0)^2)
  K_mm <- commutation_matrix(nrow(B_0), nrow(B_0))
  E_m <- elimination_matrix(nrow(B_0))
  I_n <- Matrix::Diagonal(nrow(S_0))
  I_nn <- Matrix::Diagonal(nrow(S_0)^2)
  K_nn <- commutation_matrix(nrow(S_0), nrow(S_0))
  E_n <- elimination_matrix(nrow(S_0))
  K_XTX_Sigma <- commutation_matrix(ncol(XTX), nrow(invSigma))
  I_Sigma <- Matrix::Diagonal(ncol(invSigma))
  I_Sigma2 <- Matrix::Diagonal(ncol(invSigma)^2)
  I_XTX <- Matrix::Diagonal(nrow(XTX))
  I_XTX2 <- Matrix::Diagonal(nrow(XTX)^2)
  sks_xtx <- I_Sigma %x% K_XTX_Sigma %x% I_XTX
  EIK_mm <- E_m %*% (I_mm + K_mm)
  EIK_nn <- E_n %*% (I_nn + K_nn)
  tinvB_x_invB <- t(invB) %x% invB

  pb <- txtProgressBar(1, num_steps, style = 3)
  setTxtProgressBar(pb, 0)
  for (i in 1:total) {
    # now <- Sys.time()
    # Update beta - part 1
    invSigma_x_XTX <- kronecker(invSigma, XTX)
    Vg <- solve(invB + invSigma_x_XTX)

    d_invSigma_x_XTX <- d_kronecker(invSigma, d_invSigma, XTX, d_XTX, sks_xtx,
                                    I_Sigma2, I_XTX2)
    d_invB <- d_inv(invB, d_B0, tinvB_x_invB)
    d_invB_p_invSigma_x_XTX <- d_sum(d_invB, d_invSigma_x_XTX)
    inv_of_invB_p_invSigma_x_XTX <- Vg
    d_Vg <- d_inv(inv_of_invB_p_invSigma_x_XTX, d_invB_p_invSigma_x_XTX)

    # Update beta - part 2
    invVb_sum <- invVb + as.numeric(XTY %*% invSigma)
    bg <- Vg %*% invVb_sum

    d_vec_XTY_invSigma <- d_product(XTY, d_XTY, invSigma, d_invSigma)
    d_invVb <- d_product(invB, d_invB, b_0, d_b0)
    d_invVb_sum <- d_sum(d_invVb, d_vec_XTY_invSigma)
    d_bg <- d_product(Vg, d_Vg, invVb_sum, d_invVb_sum)

    # Update beta - part 3
    chol_Vg <- t(chol(Vg))
    beta_g <- bg + chol_Vg %*% Z[,i]

    d_chol_Vg <- d_chol(chol_Vg, d_Vg, I_mm, K_mm, I_m, E_m, EIK_mm)
    d_beta_g <- d_sum(d_bg, d_product(chol_Vg, d_chol_Vg, Z[,i], d_Z_i))

    # Update Sigma
    B <- matrix(beta_g, K, n)
    tB <- t(B)
    tB_XTY <- tB %*% XTY
    XTX_B <- XTX %*% B
    S <- S_0 + YTY - tB_XTY - t(tB_XTY) + tB %*% XTX_B
    inv_S <- solve(S)
    L <- t(chol(inv_S))
    A <- diag(C[, i]) + get_lower_tri(Z2[,,i])
    LA <- L %*% A
    invSigma <- LA %*% t(LA)
    # inv_LA <- forwardsolve(LA, I_n)
    # Sigma <- t(inv_LA) %*% inv_LA
    Sigma <- solve(invSigma)

    d_B <- d_beta_g
    d_tB <- d_transpose(B, d_B)
    d_tB_XTY <- d_product(tB, d_tB, XTY, d_XTY)
    d_S <- d_S0 %>% d_sum(d_YTY) %>% d_minus(d_tB_XTY) %>%
      d_minus(d_transpose(tB_XTY, d_tB_XTY)) %>%
      d_sum(d_product(tB, d_tB, XTX_B, d_product(XTX, d_XTX, B, d_B)))
    d_A <- init_deriv(length(A))
    d_A$d_nu0 <- matrix(as.numeric(diag(
      purrr::map_dbl(1:n, function(.x) {
        d_ii <- C[.x, i]
        delta_i <- d_ii ^ 2
        delta_bar_i <- delta_i / 2
        (2 * delta_bar_i)^(-0.5) *
          0.5 * d_Gamma(delta_bar_i, 0.5 * (v_1 - .x + 1))
      })
    )))
    d_inv_S <- d_inv(inv_S, d_S)
    d_L <- d_chol(L, d_inv_S, I_nn, K_nn, I_n, E_n, EIK_nn)
    d_LA <- d_product(L, d_L, A, d_A)
    d_invSigma <- d_XXT(LA, d_LA)
    d_Sigma_on_d_invSigma <- t(Sigma) %x% Sigma

    # Keep track
    res_beta[i, ] <- beta_g
    res_sigma[i, ] <- as.numeric(Sigma)
    runs_d_beta[[i]] <- d_beta_g
    runs_d_Sigma[[i]] <- apply_chain(
      . %>% {d_Sigma_on_d_invSigma %*% d_invSigma[[.]]},
      names(d_invSigma)
    )
    setTxtProgressBar(pb, i)
    # print(sprintf("\n Estimate time to completion: %f minutes\n",
            # round(as.numeric(Sys.time() - now, units = "mins") * (total - i), 3)))
  }
  # Return
  list(beta = res_beta,
       Sigma = res_sigma,
       d_beta = runs_d_beta %>% collect_and_reshape(),
       d_Sigma = runs_d_Sigma %>% collect_and_reshape())
}
