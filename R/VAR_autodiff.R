#' Sensitivity analysis for Vector-Auto-Regressive model with independent Normal-Wishart priors.
#' @param data0 Matrix; each row is one observation, each column is one measurement / predictor.
#' @param lag Integer; the lag of the time series model.
#' @param b_0 A numeric vector; the mean for the multivariate normal prior.
#' @param B_0 A numeric matrix; the covariance for the multivariate normal prior.
#' @param v_0 scalar; degree of freedom for the inverse-Wishart prior.
#' @param S_0 matrix; the scale matrix for the inverse-Wishart prior.
#' @param init_Sigma (Optional) matrix; the starting value of the noise inverse-covariance.
#' @param num_steps integer; number of MCMC steps.
#' @param burn_ins integer; number of burn-ins.
#' @export
VAR_AD <- function(data0, lag, b_0, B_0, v_0, S_0, init_Sigma,
                     num_steps = 3e3, burn_ins = 1e3) {
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
  invSigma <- init_Sigma

  XTX <- t(X) %*% X
  XTY <- t(X) %*% Y
  YTY <- t(Y) %*% Y

  total <- num_steps + burn_ins
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
  d_XTX <- init_deriv(length(XTX))
  d_XTY <- init_deriv(length(XTY))
  d_YTY <- init_deriv(length(YTY))
  d_Z_i <- init_deriv(length(Z[,1]))
  runs_d_beta <- vector("list", total)
  runs_d_Sigma <- vector("list", total)
  I_m <- Matrix::Diagonal(nrow(B_0))
  I_mm <- Matrix::Diagonal(nrow(B_0)^2)
  K_mm <- as(matrixcalc::commutation.matrix(nrow(B_0), nrow(B_0)), "dgCMatrix")
  E_m <- as(matrixcalc::elimination.matrix(nrow(B_0)), "dgCMatrix")
  I_n <- Matrix::Diagonal(nrow(S_0))
  I_nn <- Matrix::Diagonal(nrow(S_0)^2)
  K_nn <- as(matrixcalc::commutation.matrix(nrow(S_0), nrow(S_0)), "dgCMatrix")
  E_n <- as(matrixcalc::elimination.matrix(nrow(S_0)), "dgCMatrix")
  K_XTX_Sigma <- as(matrixcalc::commutation.matrix(ncol(XTX), nrow(invSigma)), "dgCMatrix")
  I_Sigma <- Matrix::Diagonal(ncol(invSigma))
  I_XTX <- Matrix::Diagonal(nrow(XTX))

  pb <- txtProgressBar(1, num_steps, style = 3)
  for (i in 1:total) {
    # Update beta - part 1
    invSigma_x_XTX <- kronecker(invSigma, XTX)
    Vg <- solve(invB + invSigma_x_XTX)

    d_invSigma_x_XTX <- d_kronecker(invSigma, d_invSigma, XTX, d_XTX, I_Sigma, K_XTX_Sigma, I_XTX)
    d_invB <- d_inv(invB, d_B0)
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

    d_chol_Vg <- d_chol(chol_Vg, d_Vg, I_mm, K_mm, I_m, E_m)
    d_beta_g <- d_sum(d_bg, d_product(chol_Vg, d_chol_Vg, Z[,i], d_Z_i))

    # Update Sigma
    B <- matrix(beta_g, K, n)
    tB <- t(B)
    tB_XTY <- tB %*% XTY
    S <- S_0 + YTY - tB_XTY - t(tB_XTY) + tB %*% XTX %*% B
    inv_S <- solve(S)
    L <- t(chol(inv_S))
    A <- diag(C[, i]) + get_lower_tri(Z2[,,i])
    LA <- L %*% A
    invSigma <- LA %*% t(LA)
    # inv_LA <- forwardsolve(LA, I_n)
    # Sigma <- t(inv_LA) %*% inv_LA
    Sigma <- solve(invSigma)

    d_B <- d_beta_g
    d_tB_XTY <- d_product(tB, d_transpose(B, d_B), XTY, d_XTY)
    d_S <- d_S0 %>% d_sum(d_YTY) %>% d_minus(d_tB_XTY) %>%
      d_minus(d_transpose(tB_XTY, d_tB_XTY)) %>%
      d_sum(d_product(
        tB, d_transpose(B, d_B),
        XTX %*% B, d_product(XTX, d_XTX, B, d_B)
      ))
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
    d_L <- d_chol(L, d_inv_S, I_nn, K_nn, I_n, E_n)
    d_LA <- d_product(L, d_L, A, d_A)
    d_invSigma <- d_XXT(LA, d_LA)
    d_Sigma_on_d_invSigma <- t(Sigma) %x% Sigma

    # Keep track
    res_beta[i, ] <- beta_g
    res_sigma[i, ] <- as.numeric(Sigma)
    runs_d_beta[[i]] <- d_beta_g
    runs_d_Sigma[[i]] <- apply_chain(. %>% {d_Sigma_on_d_invSigma %*% d_invSigma[[.]]})
    setTxtProgressBar(pb, i)
  }
  # Return
  list(beta = res_beta %>% tail(num_steps),
       Sigma = res_sigma %>% tail(num_steps),
       d_beta = runs_d_beta %>% tail(num_steps) %>% collect(),
       d_Sigma = runs_d_Sigma %>% tail(num_steps) %>% collect())
}


apply_chain <- function(expr_fun) {
  hyperparameter <- c("d_b0", "d_B0", "d_nu0", "d_S0", "d_Sigma0")
  hyperparameter %>%
    purrr::map(expr_fun) %>%
    purrr::set_names(hyperparameter)
}


d_chol <- function(L, dA, I_nn, K_nn, I_n, E_n) {
  # LL^T = A
  n <- nrow(L)
  if (missing(I_n)) I_n <- diag(n)
  if (missing(I_nn)) I_nn <- diag(n^2)
  if (missing(K_nn)) K_nn <- matrixcalc::commutation.matrix(n, n)
  if (missing(E_n)) E_n <- matrixcalc::elimination.matrix(n)
  D_n <- Matrix::t(E_n)
  fac_1 <- E_n %*% (I_nn + K_nn)
  fac_2 <- D_n %*% solve(fac_1 %*% kronecker_sp_3_cpp(L, as.matrix(D_n))) %*% E_n
  # apply_chain(. %>% {fac_2 %*% dA[[.]]})
  apply_chain(. %>% {eigenMapMatMult(as.matrix(fac_2), dA[[.]])})
}
d_transpose <- function(X, dX, K_nq) {
  n <- nrow(X)
  q <- ncol(X)
  if (missing(K_nq)) K_nq <- matrixcalc::commutation.matrix(n, q)
  apply_chain(. %>% {K_nq %*% dX[[.]]})
}
d_XXT <- function(X, dX, I_nn, K_nn, I_n) {
  n <- nrow(X)
  if (missing(I_n)) I_n <- diag(n)
  if (missing(I_nn)) I_nn <- diag(n^2)
  if (missing(K_nn)) K_nn <- matrixcalc::commutation.matrix(n, n)
  apply_chain(. %>% {(I_nn + K_nn) %*% (X %x% I_n) %*% dX[[.]]})
}
d_kronecker <- function(A, dA, B, dB, I_n, K_qm, I_m) {
  # apply_chain(. %>% {A %x% dB[[.]] + dA[[.]] %x% B})
  m <- nrow(A)
  n <- ncol(A)
  p <- nrow(B)
  q <- ncol(B)
  if (missing(I_n)) I_n <- diag(n)
  if (missing(K_qm)) K_qm <- matrixcalc::commutation.matrix(q, m)
  if (missing(I_m)) I_m <- diag(m)
  fac_1 <- (I_n %x% K_qm %x% I_m)
  apply_chain(. %>% {fac_1 %*%
      (as.numeric(A) %x% dB[[.]] + dA[[.]] %x% as.numeric(B))})
}
d_sum <- function(dA, dB) {
  apply_chain(. %>% {dA[[.]] + dB[[.]]})
}
d_product <- function(A, dA, B, dB, I_a, I_b) {
  if (missing(I_a)) I_a <- diag(ifelse(is.vector(A), length(A), nrow(A)))
  if (missing(I_b)) I_b <- diag(ifelse(is.vector(B), 1, ncol(B)))
  apply_chain(. %>% {(I_b %x% A) %*% dB[[.]] + (t(B) %x% I_a) %*% dA[[.]]})
}
d_minus <- function(dA, dB) {
  apply_chain(. %>% {dA[[.]] - dB[[.]]})
}
d_inv <- function(inv_A, dA) {
  fac_1 <- t(inv_A) %x% inv_A
  # apply_chain(. %>% {fac_1 %*% dA[[.]]})
  apply_chain(. %>% {eigenMapMatMult(fac_1, as.matrix(dA[[.]]))})
}
d_Gamma <- function(g, alpha) {
  f <- function(t) { log(t) * dgamma(t, alpha, 1) }
  num_1 <- integrate(f, 0, g)$value
  num_2 <- digamma(alpha) * pgamma(g, alpha, 1)
  - (num_1 - num_2) / dgamma(g, alpha, 1)
}
