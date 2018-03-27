#' Sensitivity analysis for Vector-Auto-Regressive model with normal-Wishart priors.
#' @param data0 A numeric matrix; a multivariate time series. The number
#' of columns is equal to the number of time points, the number of rows is
#' equal to the dimension of the observation. One can get such data from
#' the first component of the output from 'simulate_VAR'.
#' @param p integer; lag of the time series model.
#' @param b_0 A numeric vector; the mean for the multivariate normal prior.
#' @param V A numeric matrix; the covariance for the multivariate normal prior.
#' @param v_0 scalar; degree of freedom for the inverse-Wishart prior.
#' @param S_0 matrix; the scale matrix for the inverse-Wishart prior.
#' @param init_Sigma (Optional) matrix; the starting value of the noise covariance.
#' @param num_steps integer; number of MCMC steps.
#' @param burn_ins integer; number of burn-ins.
#' @export
VAR_AD <- function(data0, p, b_0, V, v_0, S_0, init_Sigma,
                   num_steps = 1e4, burn_ins = 1e3) {
  if (missing(init_Sigma))
    init_Sigma <- MCMCpack::riwish(v_0, S_0)

  # Initialisation
  T0 <- ncol(data0) - p    # number of data
  n <- nrow(data0)         # dimension of data
  len_beta <- length(b_0)  # dimension of parameter q
  y <- matrixcalc::vec(data0[, seq(T0) + p])
  X <- VAR_model_matrix(data0, p)
  inv_V <- solve(V)
  Sigma_g <- init_Sigma
  inv_V_times_b_0 <- inv_V %*% b_0
  v_1 <- v_0 + T0
  keep <- num_steps - burn_ins
  runs_param <- vector("list", num_steps)
  runs_d_beta <- vector("list", num_steps)
  runs_d_Sigma <- vector("list", num_steps)

  # Autodiff helper variables and functions
  I_T <- Matrix::Diagonal(T0)
  I_n <- Matrix::Diagonal(n)
  I_q <- Matrix::Diagonal(len_beta)
  I_nn <- Matrix::Diagonal(n^2)
  I_qq <- Matrix::Diagonal(len_beta^2)
  K_nn <- as(matrixcalc::commutation.matrix(n, n), "dgCMatrix")
  K_qq <- as(matrixcalc::commutation.matrix(len_beta), "dgCMatrix")
  K_nT <- as(matrixcalc::commutation.matrix(n, T0), "dgCMatrix")
  elimL <- as(matrixcalc::elimination.matrix(len_beta), "dgCMatrix")
  E <- as(matrixcalc::elimination.matrix(n), "dgCMatrix")
  D <- Matrix::t(E)

  vec_I_T <- as(matrix(as.numeric(I_T)), "dgCMatrix")
  tX <- t(X)
  XT_otimes_XT <- as(tX %x% tX, "dgCMatrix")
  y_T_otimes_X_T <- as(t(y) %x% tX, "dgCMatrix")
  IT_KnT_In <- I_T %x% K_nT %x% I_n
  n_ivt_otimes_iv <- neg_tx_otimes_x(inv_V)
  b_0_T_otimes_I_q <- as(t(b_0) %x% I_q, "dgCMatrix")
  d_bg_fac_2 <- b_0_T_otimes_I_q %*% n_ivt_otimes_iv
  I_nn_p_K_nn <- I_nn + K_nn
  I_qq_p_K_qq <- I_qq + K_qq
  elimL_x_I_qq_p_K_qq <- elimL %*% (I_qq_p_K_qq)
  E_I_nn_p_K_nn <- E %*% I_nn_p_K_nn

  hyperparameter <- c("d_b0", "d_V", "d_nu0", "d_S0", "d_Sigma0")
  d_Sigma <- init_VAR_differential(n^2, len_beta, n)
  d_Sigma$d_Sigma0 <- I_nn

  apply_chain <- function(expr_fun) {
    purrr::set_names(
      purrr::map(hyperparameter, ~expr_fun(.x)),
      hyperparameter
    )
  }
  deriv_Ag <- function(inv_Sigma_g, d_Sigma) {
    fac_2 <- neg_tx_otimes_x(inv_Sigma_g)
    expr <- . %>% {IT_KnT_In %*% (vec_I_T %x% (fac_2 %*% d_Sigma[[.]]))}
    d_Ag <- apply_chain(expr)
    d_Ag
  }
  deriv_Kg <- function(d_Ag) {
    expr <- . %>% {XT_otimes_XT %*% d_Ag[[.]]}
    d_Kg <- apply_chain(expr)
    d_Kg$d_V <- d_Kg$d_V + n_ivt_otimes_iv
    d_Kg
  }
  deriv_bg <- function(inv_Sigma_g, inv_K_g, d_Kg, d_Ag) {
    # A_g <- (I_T %x% inv_Sigma_g)
    # fac_1 <- (Matrix::t(inv_V_times_b_0 + tX %*% A_g %*% y) %x% I_q) %*% neg_tx_otimes_x(inv_K_g)
    fac_1 <- kronecker_sp_3_cpp(
      t(as.numeric(inv_V_times_b_0 + kronecker_sp_1(tX, inv_Sigma_g) %*% y)),
      neg_tx_otimes_x(inv_K_g)
    )
    fac_2 <- d_bg_fac_2
    fac_3 <- y_T_otimes_X_T
    expr <- . %>% {fac_1 %*% d_Kg[[.]] + inv_K_g %*% (fac_3 %*% d_Ag[[.]])}
    d_bg <- apply_chain(expr)
    d_bg$d_b0 <- d_bg$d_b0 + inv_K_g %*% inv_V
    d_bg$d_V <- d_bg$d_V + inv_K_g %*% fac_2
    d_bg
  }
  deriv_beta <- function(inv_Sigma_g, d_Sigma, inv_K_g, L, z) {
    # fac_1 <- Matrix::t(elimL) %*%
    #   solve(elimL %*% (I_qq_p_K_qq) %*% (L %x% I_q) %*% Matrix::t(elimL)) %*%
    #   elimL
    fac_1 <- Matrix::t(elimL) %*%
      solve(
        elimL_x_I_qq_p_K_qq %*%
          kronecker_sp_3_cpp(as.matrix(L), as.matrix(Matrix::t(elimL)))
      ) %*%
      elimL
    # fac_2 <- (t(z) %x% I_q) %*% fac_1 %*% neg_tx_otimes_x(inv_K_g)
    fac_2 <- kronecker_sp_3_cpp(t(z), as.matrix(fac_1)) %*% neg_tx_otimes_x(inv_K_g)
    d_Ag <- deriv_Ag(inv_Sigma_g, d_Sigma)
    d_Kg <- deriv_Kg(d_Ag)
    d_bg <- deriv_bg(inv_Sigma_g, inv_K_g, d_Kg, d_Ag)

    expr <- . %>% {d_bg[[.]] + fac_2 %*% d_Kg[[.]]}
    d_beta <- apply_chain(expr)
    d_beta
  }

  deriv_G <- function(G, alpha, update_var) {
    d_G <- init_VAR_differential(1, len_beta, n)
    f <- function(t) { log(t) * dgamma(t, alpha, 1) }
    num_1 <- integrate(f, 0, G)$value
    num_2 <- digamma(alpha) * pgamma(G, alpha, 1)
    # 0.5 is the factor associated with alpha at function call
    d_G[[update_var]] <- - 0.5 * (num_1 - num_2) / dgamma(G, alpha, 1)
    d_G
  }
  deriv_inv_d_ii <- function(delta_i, dg) {
    d_G <- deriv_G(delta_i, dg / 2, "d_nu0")
    d_G$d_nu0 <- - (2 * delta_i)^(-1.5) * d_G$d_nu0
    d_G
  }
  deriv_B <- function() {
    dg <- v_1 - seq(n) + 1
    delta_i <- rgamma(n, dg / 2, 1)
    list_d_inv_d_ii <- purrr::map2(delta_i, dg, .f = ~deriv_inv_d_ii(.x, .y))
    d_B <- init_VAR_differential(n^2, len_beta, n)
    d_B$d_nu0[seq(1, n^2, n+1), ] <- purrr::map_dbl(list_d_inv_d_ii, ~.$d_nu0)
    d_B
  }
  find_sum_term <- function(beta_g) {
    m_vec <- y - X %*% beta_g
    res <- 0
    for (t in 1:T0) {
      ls <- (1 + n*(t-1)) : (n*t)
      res <- res + kronecker_sp_2(matrix(m_vec[ls]), -X[ls,]) +
        kronecker_sp_3_cpp(matrix(m_vec[ls]), -X[ls,])
    }
    res
  }
  deriv_Sigma <- function(beta_g, d_beta, R, B) {
    d_B <- deriv_B()
    sum_term <- find_sum_term(beta_g)
    RB <- R %*% B
    # fac_1 <- I_nn_p_K_nn %*% ((RB %*% t(B)) %x% I_n) %*%
    #   (D %*% solve(E %*% I_nn_p_K_nn %*% (R %x% I_n) %*% D) %*% E) %*%
    #   sum_term
    fac_1_core <- I_nn_p_K_nn %*%
      kronecker_sp_3_cpp(
        as.matrix(RB %*% t(B)),
        as.matrix((D %*% solve(E_I_nn_p_K_nn %*%
                       kronecker_sp_3_cpp(as.matrix(R), as.matrix(D))) %*% E))
      )
    fac_1 <- fac_1_core %*% sum_term
    fac_2 <- I_nn_p_K_nn %*% (RB %x% R)
    expr <- . %>% {fac_1 %*% d_beta[[.]] + fac_2 %*% d_B[[.]]}
    d_Sigma <- apply_chain(expr)
    fac_1b <- fac_1_core %*% (I_nn + sum_term %*% d_beta$d_S0)
    d_Sigma$d_S0 <- fac_1b  + fac_2 %*% d_B$d_S0
    d_Sigma
  }
  Bartlett_A <- function(v_1, n) {
    m0 <- matrix(0, n, n)
    m0[lower.tri(m0)] <- rnorm(n * (n-1) / 2)
    diag(m0) <- rchisq(n, df = v_1 - 1:n + 1)
    m0
  }

  pb <- txtProgressBar(1, num_steps, style = 3)
  for (i in 1:num_steps) {
    # Update beta
    inv_Sigma_g <- solve(Sigma_g)
    tX_fac <- kronecker_sp_1(tX, inv_Sigma_g)
    K_g <- inv_V + tX_fac %*% X
    inv_K_g <- solve(K_g)
    b_g <- inv_K_g %*% (inv_V_times_b_0 + tX_fac %*% y)
    L <- t(chol(inv_K_g))
    z <- rnorm(len_beta)
    beta_g <- b_g + L %*% z
    # Autodiff beta
    d_beta <- deriv_beta(inv_Sigma_g, d_Sigma, inv_K_g, L, z)

    # Update Sigma
    S_g <- S_0 + residuals_cov(y, X, beta_g, n)
    R <- t(chol(S_g))
    A <- Bartlett_A(v_1, nrow(S_g))
    inv_A <- forwardsolve(A, diag(n))
    Sigma_g <- R %*% t(inv_A) %*% inv_A %*% t(R)
    # Autodiff Sigma
    d_Sigma <- deriv_Sigma(beta_g, d_beta, R, t(inv_A))

    # Keep track
    runs_param[[i]] <- list(beta = as.numeric(beta_g), Sigma = Sigma_g)
    runs_d_beta[[i]] <- d_beta
    runs_d_Sigma[[i]] <- d_Sigma
    setTxtProgressBar(pb, i)
  }

  # Tidy format
  tidy_format <- . %>% tail(keep) %>% do.call(rbind, .)
  list(
    Sigma = runs_param %>% purrr::map(~as.numeric(.x$Sigma)) %>% tidy_format(),
    beta = runs_param %>% purrr::map(~Matrix::t(.x$beta)) %>% tidy_format(),
    d_Sigma = runs_d_Sigma %>% tail(keep) %>% tidy_VAR_differential(),
    d_beta = runs_d_beta %>% tail(keep) %>% tidy_VAR_differential()
  )
}


init_VAR_differential <- function(len0, q, n) {
  # len0 is the number of entries in the variable of interest.
  list(
    d_b0 = zeros(len0, q),
    d_V = zeros(len0, q^2),
    d_nu0 = zeros(len0, 1),
    d_S0 = zeros(len0, n^2),
    d_Sigma0 = zeros(len0, n^2)
  )
}


tidy_VAR_differential <- function(dlist0) {
  extract_rbind <- function(attr0) {
    dlist0 %>%
      purrr::map(~as.numeric(.x[[attr0]])) %>%
      do.call(rbind, .)
  }
  list(
    d_b0 = extract_rbind("d_b0"),
    d_V = extract_rbind("d_V"),
    d_nu0 = extract_rbind("d_nu0"),
    d_S0 = extract_rbind("d_S0"),
    d_Sigma0 = extract_rbind("d_Sigma0")
  )
}
