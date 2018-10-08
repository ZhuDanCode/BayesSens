#' Sensitivity analysis for 2-equation SUR (seemingly unrelated regression) model with multivariate normal
#' prior for the means and inverse gamma for the variance.
#' @param Xy A numeric matrix; the first set of covariates.
#' @param y A numeric vector; the response variable for the first set of covariates.
#' @param b_0 A numeric vector; the mean for the first multivariate normal prior.
#' @param B_0 A numeric matrix; the covariance for the first multivariate normal prior.
#' @param Xs A numeric matrix; the second set of covariates.
#' @param s A numeric vector; the response variable for the second set of covariates.
#' @param g_0 A numeric vector; the mean for the second multivariate normal prior.
#' @param G_0 A numeric matrix; the covariance for the second multivariate normal prior.
#' @param v_0 A positive number; the degrees of freedom of the Wishart prior.
#' @param R_0 A 2x2 symmetric matrix; the scale matrix of the Wishart prior.
#' @param init_gamma (Optional) A numeric vector; the starting value of the gamma parameter.
#' @param init_Sigma (Optional) A 2x2 symmetric matrix; the starting value of the sigma parameter.
#' @param num_steps Integer; number of MCMC steps.
#' @examples
#' \dontrun{
#' n <- 1000
#' p <- 2
#' k <- 3
#' data0 <- SUR2_data(n, p, k,
#'                       Sigma = matrix(c(1, 0.2, 0.2, 1), 2, 2),
#'                       intercept_1 = TRUE, intercept_2 = TRUE)
#' res <- SUR2_AD(Xy = data0$Xy, y = data0$y,
#'                   b_0 = numeric(p + 2), B_0 = diag(p + 2),
#'                   Xs = data0$Xs, s = data0$s,
#'                   g_0 = numeric(k + 1), G_0 = diag(k + 1),
#'                   v_0 = 5, R_0 = diag(2),
#'                   num_steps = 1e3)
#' }
#' @export
#' @references Jacobi, Liana and Joshi, Mark S. and Zhu, Dan, Automated Sensitivity Analysis for Bayesian Inference via Markov Chain Monte Carlo: Applications to Gibbs Sampling (February 9, 2018). Available at SSRN: https://ssrn.com/abstract=2984054 or http://dx.doi.org/10.2139/ssrn.2984054
SUR2_AD <- function(Xy, y, b_0, B_0, Xs, s, g_0, G_0, v_0, R_0,
                       init_gamma, init_Sigma, num_steps = 1e4) {
  if (missing(init_gamma))
    init_gamma <- numeric(length(g_0))
  if (missing(init_Sigma))
    init_Sigma <- diag(2)

  # Initialisation
  n <- length(y)
  v_1 <- v_0 + n
  gamma <- init_gamma
  Sigma <- solve(init_Sigma)

  res_beta <- matrix(0, num_steps,length(b_0))
  res_gamma <- matrix(0, num_steps,length(g_0))
  res_Sigma <- matrix(0, num_steps,4)
  runs_param <- vector("list", num_steps)
  runs_d_beta <- vector("list", num_steps)
  runs_d_gamma <- vector("list", num_steps)
  runs_d_Sigma <- vector("list", num_steps)

  len_b0 <- length(b_0)
  len_g0 <- length(g_0)

  #Precomputations
  invB_0 <- solve(B_0)
  invB_0_b_0 <- invB_0 %*% b_0
  invG_0 <- solve(G_0)
  invG_0_g_0 <- invG_0 %*% g_0
  # invR_0 <- solve(R_0)
  tinvB_0_x_invB_0 <- t(invB_0) %x% invB_0
  tinvG_0_x_invG_0 <- t(invG_0) %x% invG_0

  v <- cbind(y, s)
  X <- cbind(Xy, Xs)

  XyTXy <- crossprod(Xy, Xy)
  XyTXs <- crossprod(Xy, Xs)
  XsTXy <- crossprod(Xs, Xy)
  XsTXs <- crossprod(Xs, Xs)
  XyTy <- crossprod(Xy, y)
  XyTs <- crossprod(Xy, s)
  XsTy <- crossprod(Xs, y)
  XsTs <- crossprod(Xs, s)
  vTv <- crossprod(v, v)
  vTX <- crossprod(v, X)

  #pre-simluation
  Zy <- matrix(rnorm(len_b0 * num_steps), len_b0, num_steps)
  Zs <- matrix(rnorm(len_g0 * num_steps), len_g0, num_steps)
  Z2 <- array(rnorm(2 * 2 * num_steps), dim = c(2, 2, num_steps)) #2 because there are 2 equations
  C <- matrix(0, 2, num_steps) #2 because there are 2 equations
  for (i in 1:2) {
    C[i,] <- rchisq(num_steps, v_1 - i + 1)
  }
  C <- sqrt(C)

  #Cholesky
  I_p <- Matrix::Diagonal(len_b0)
  I_pp <- Matrix::Diagonal(len_b0^2)
  K_pp <- commutation_matrix(len_b0, len_b0)
  E_p <- elimination_matrix(len_b0)
  EIK_pp <- E_p %*% (I_pp + K_pp)
  I_k <- Matrix::Diagonal(len_g0)
  I_kk <- Matrix::Diagonal(len_g0^2)
  K_kk <- commutation_matrix(len_g0, len_g0)
  E_k <- elimination_matrix(len_g0)
  EIK_kk <- E_k %*% (I_kk + K_kk)

  I_2 <- diag(2)
  I_22 <- diag(4)
  K_22 <- commutation_matrix(2, 2)
  E_2 <- elimination_matrix(2)
  EIK_22 <- E_2 %*% (I_22 + K_22)

  #Autodiff
  init_deriv <- . %>% {list(
    d_gamma0 = matrix(0, ., length(init_gamma)),
    d_Sigma0 = matrix(0, ., length(init_Sigma)),
    d_b0 = matrix(0, ., length(b_0)),
    d_B0 = matrix(0, ., length(B_0)),
    d_g0 = matrix(0, ., length(g_0)),
    d_G0 = matrix(0, ., length(G_0)),
    d_v0 = matrix(0, ., length(v_0)),
    d_R0 = matrix(0, ., length(R_0))
  )}

  d_b0 <- init_deriv(length(b_0))
  diag(d_b0$d_b0) <- 1
  d_B0 <- init_deriv(length(B_0))
  d_B0$d_B0 <- differential_matrix(nrow(B_0))
  d_g0 <- init_deriv(length(g_0))
  diag(d_g0$d_g0) <- 1
  d_G0 <- init_deriv(length(G_0))
  d_G0$d_G0 <- differential_matrix(nrow(G_0))
  d_R0 <- init_deriv(length(R_0))
  d_R0$d_R0 <- differential_matrix(nrow(R_0))
  d_gamma <- init_deriv(length(init_gamma))
  diag(d_gamma$d_gamma0) <- 1
  d_Sigma <- init_deriv(length(init_Sigma))
  d_Sigma$d_Sigma0 <- differential_matrix(nrow(init_Sigma))

  d_invB_0 <- d_inv(invB_0, d_B0, - tinvB_0_x_invB_0)
  d_invB_0_b_0 <- d_product(invB_0, d_invB_0, b_0, d_b0)
  d_invG_0 <- d_inv(invG_0, d_G0, - tinvG_0_x_invG_0)
  d_invG_0_g_0 <- d_product(invG_0, d_invG_0, g_0, d_g0)
  invR_0 <- solve(R_0)
  d_invR_0 <- d_inv(invR_0, d_R0)

  d_XyTXy <- 0
  d_XyTXs <- 0
  d_XsTXy <- 0
  d_XsTXs <- 0
  d_XyTy <- 0
  d_XyTs <- 0
  d_XsTy <- 0
  d_XsTs <- 0

  d_vTv <- 0
  d_vTX <- 0
  d_X <- 0
  d_v <- 0

  d_Zy_i <- 0
  d_Zs_i <- 0
  d_Z2 <- 0

  pb <- txtProgressBar(1, num_steps, style = 3)
  setTxtProgressBar(pb, 0)
  for (i in 1:num_steps) {
    #loop
    sigma11 <- Sigma[1]
    sigma12 <- Sigma[2]
    sigma21 <- Sigma[3]
    sigma22 <- Sigma[4]
    invw11 <- 1 / (sigma11 - sigma12 * sigma21 / sigma22)
    invw22 <- 1 / (sigma22 - sigma12 * sigma21 / sigma11)
    d_sigma11 <- lapply(d_Sigma, function(x) t(matrix(x[1, ])))
    d_sigma12 <- lapply(d_Sigma, function(x) t(matrix(x[2, ])))
    d_sigma21 <- lapply(d_Sigma, function(x) t(matrix(x[3, ])))
    d_sigma22 <- lapply(d_Sigma, function(x) t(matrix(x[4, ])))
    d_sigma12_sigma21 <- d_constant_multiply_constant(
      sigma12, d_sigma12,
      sigma21, d_sigma21
    )
    d_sigma12_sigma21_on_sigma22 <- d_constant_divide_constant(
      sigma12 * sigma21, d_sigma12_sigma21,
      sigma22, d_sigma22
    )
    d_sigma12_sigma21_on_sigma11 <- d_constant_divide_constant(
      sigma12 * sigma21, d_sigma12_sigma21,
      sigma11, d_sigma11
    )
    d_invw11 <- d_constant_inv(
      1 / (sigma11 - sigma12 * sigma21 / sigma22),
      d_minus(d_sigma11, d_sigma12_sigma21_on_sigma22)
    )
    d_invw22 <- d_constant_inv(
      1 / (sigma22 - sigma12 * sigma21 / sigma11),
      d_minus(d_sigma22, d_sigma12_sigma21_on_sigma11)
    )

    #Update beta
    #Update Bg
    invw11_XyTXy <- invw11 * XyTXy
    Bg <- solve(invB_0 + invw11_XyTXy)
    #Update d_Bg
    d_invw11_XyTXy <- d_constant_multiply_matrix(invw11, d_invw11, XyTXy, d_XyTXy)
    d_Bg <- d_inv(Bg, d_sum(d_invB_0, d_invw11_XyTXy))

    #Update bg
    b2 <- XyTs - XyTXs %*% gamma
    b4 <- XyTy - sigma12 * sigma21 / sigma22 * b2
    b7 <- invB_0_b_0 + invw11 * b4
    bg <- Bg %*% b7
    #Update d_bg
    d_b2 <- d_minus(d_XyTs, d_product(XyTXs, d_XyTXs, gamma, d_gamma))
    d_b4 <- d_minus(
      d_XyTy,
      d_constant_multiply_matrix(
        sigma12 * sigma21 / sigma22, d_sigma12_sigma21_on_sigma22,
        b2, d_b2
      )
    )
    d_b7 <- d_sum(
      d_invB_0_b_0,
      d_constant_multiply_matrix(invw11, d_invw11, b4, d_b4)
    )
    d_bg <-  d_product(Bg, d_Bg, b7, d_b7)

    #Update betag
    chol_Bg <- t(chol(Bg))
    beta <- bg + chol_Bg %*% Zy[, i]
    #Update d_betag
    d_chol_Bg <- d_chol(chol_Bg, d_Bg, I_pp, K_pp, I_p, E_p, EIK_pp)
    d_beta <- d_sum(d_bg, d_product(chol_Bg, d_chol_Bg, Zy[, i], d_Zy_i))

    #Update gamma
    #Update Gg
    invw22_XsTXs <- invw22 * XsTXs
    Gg <- solve(invG_0 + invw22_XsTXs)
    #Update d_Gg
    d_invw22_XsTXs <- d_constant_multiply_matrix(
      invw22, d_invw22, XsTXs, d_XsTXs
    )
    d_Gg <- d_inv(Gg, d_sum(d_invG_0, d_invw22_XsTXs))

    #Update gg
    g2 <- XsTy - XsTXy %*% beta
    g4 <- XsTs - sigma12 * sigma21 / sigma11 * g2
    g7 <- invG_0_g_0 + invw22 * g4
    gg <- Gg %*% g7
    #Update d_gg
    d_g2 <- d_minus(d_XsTy, d_product(XsTXy, d_XsTXy, beta, d_beta))
    d_g4 <- d_minus(d_XsTs, d_constant_multiply_matrix(
      sigma12 * sigma21 / sigma11, d_sigma12_sigma21_on_sigma11,
      g2, d_g2
    ))
    d_g7 <- d_sum(
      d_invG_0_g_0,
      d_constant_multiply_matrix(invw22, d_invw22, g4, d_g4)
    )
    d_gg <- d_product(Gg, d_Gg, g7, d_g7)

    #Update gammag
    chol_Gg <- t(chol(Gg))
    gamma <- gg + chol_Gg %*% Zs[, i]
    #Update d_gammag
    d_chol_Gg <- d_chol(chol_Gg, d_Gg, I_kk, K_kk, I_k, E_k, EIK_kk)
    d_gamma <- d_sum(d_gg, d_product(chol_Gg, d_chol_Gg, Zs[, i], d_Zs_i))

    #Update Sigma
    #Update Rg
    delta <- matrix(c(beta, rep(0, len_b0 + len_g0), gamma), ncol = 2)
    vTX_delta <- vTX %*% delta
    X_delta <- X %*% delta
    R <- invR_0 + vTv - vTX_delta - t(vTX_delta) + crossprod(X_delta)
    inv_R <- solve(R)

    #Update d_Rg
    d_delta <- add_zeros(d_beta, len_b0 + len_g0, d_gamma)
    R1 <- v - X_delta
    R2 <- t(R1) %*% R1
    d_R1 <- d_minus(d_v, d_product(X, d_X, delta, d_delta))
    d_R2 <- d_XXT(t(R1), d_transpose(R1, d_R1))
    d_R <- d_sum(d_R2, d_invR_0)
    d_inv_R <- d_inv(as.matrix(inv_R), d_R)

    #Update Sigmag
    L <- t(chol(as.matrix(inv_R)))
    A <- diag(C[, i]) + get_lower_tri(Z2[,,i])
    LA <- L %*% A
    invSigma <- LA %*% t(LA)
    Sigma <- solve(invSigma)

    #Update d_Sigmag
    d_L <- d_chol(L, d_inv_R, I_22, K_22, I_2, E_2, EIK_22)
    d_A <- init_deriv(length(A))
    d_A$d_v0 <- matrix(as.numeric(diag(
      purrr::map_dbl(1:2, function(.x) { #2 because there are 2 equations
        d_ii <- C[.x, i]
        delta_i <- d_ii ^ 2
        delta_bar_i <- delta_i / 2
        (2 * delta_bar_i)^(-0.5) *
          0.5 * d_Gamma(delta_bar_i, 0.5 * (v_1 - .x + 1))
      })
    )))
    d_invSigma <- d_XXT(LA, d_product(L, d_L, A, d_A))
    d_Sigma <- d_inv(Sigma, d_invSigma)

    # Keep track
    res_beta[i, ] <- beta
    res_gamma[i, ] <- gamma
    res_Sigma[i, ] <- as.numeric(Sigma)
    runs_d_beta[[i]] <- d_beta
    runs_d_gamma[[i]] <- d_gamma
    runs_d_Sigma[[i]] <- d_Sigma
    setTxtProgressBar(pb, i)
  }

  # Return
  return(list(beta = res_beta,
              gamma = res_gamma,
              Sigma = res_Sigma,
              d_beta = tidy_list(runs_d_beta),
              d_gamma = tidy_list(runs_d_gamma),
              d_Sigma = tidy_list(runs_d_Sigma)))
}

add_zeros <- function(a, n, b) {
  # function to stack elements of list a and b with n zeros in between
  purrr::map(
    1:length(a),
    function(x) {
      append(list(a[[x]]), as.vector(rep(0, n), mode='list')) %>%
      append(., list(b[[x]])) %>%
      do.call(rbind, .)
    }
  ) %>%
    purrr::set_names(names(a))
}

get_lower_tri <- . %>% {.[upper.tri(., T)] <- 0; .}
