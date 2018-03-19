#' Simulation data from a Vector AutoRegressive (VAR) model
#' @param m integer; number of datapoints.
#' @param n integer; dimension of the observation.
#' @param p integer; the time-series lag.
#' @param init_y_t numeric matrix (optional); the initial state of the time series.
#' The matrix must have as many columns as the lag, and the first
#' column should be the state at the first time point.
#' @param b_0 numeric vector (optional); the intercept.
#' @param B list of numeric matrices (optional); the regression coefficients.
#' @param Sigma numeric matrix (optional); the covariance matrix of the error.
#' @examples
#' simulate_VAR(30, 3, 2)
#' @export
simulate_VAR <- function(m, n, p, init_y_t, b_0, B, Sigma) {
  if (missing(b_0)) b_0 <- rnorm(n)
  if (missing(B)) B <- purrr::map(1:p, ~stationary_matrix(n, 1/p))
  if (missing(Sigma)) Sigma <- pdmatrix(n)$Sigma

  y_t <- matrix(0, nrow = n, ncol = m)
  if (missing(init_y_t)) {
    y_t[,(m-p+1):m] <- matrix(rnorm(n*p), nrow = n, ncol = p)
  } else {
    assertthat::are_equal(ncol(init_y_t), p)
    y_t[,(m-p+1):m] <- init_y_t
  }

  for (i in (m-p):1) {
    y_t[,i] <- project_ts(b_0, B, y_t[, (i+1):(i+p), drop = FALSE], Sigma)
  }

  list(data = y_t[, rev(1:m)], model = list(b_0 = b_0, B = B, Sigma = Sigma))
}

project_ts <- function(b, B, y_t, Sigma) {
  z <- b
  for (i in seq_along(B)) {
    z <- z + B[[i]] %*% y_t[,i]
  }
  z + MASS::mvrnorm(mu = numeric(length(b)), Sigma = Sigma)
}

stationary_matrix <- function(k, upper_bound = 1) {
  res <- svd(pdmatrix(k)$Sigma)
  res$u %*% diag(upper_bound * res$d / max(res$d)) %*% t(res$v)
}


# VAR util functions
#' Vectorise all the coefficients appeared in the model
#' @param model0 A VAR model; the second component of the output from 'simulate_VAR'.
#' @export
VAR_model_to_vec <- function(model0) {
  matrixcalc::vec(t(cbind(model0$b_0, do.call(cbind, model0$B))))
}

#' Convert vector back to model coefficients
#' @param vec0 Numeric vector; vectorised model coefficients.
#' @param n integer; dimension of the data.
#' @export
VAR_vec_to_model_coeff <- function(vec0, n) {
  m0 <- matrix(vec0, ncol = n, byrow = F)
  s <- seq(2, nrow(m0), n)
  list(
    b_0 = m0[1,],
    B = purrr::map2(s, s+n-1, ~t(m0[.x:.y, , drop = F]))
  )
}


#' Convert the data into a format ready for modelling
#' @param data0 Matrix; the number of columns is equal to the number of time points,
#' the number of rows is equal to the dimension of the observation. One can get
#' such data from the first component of the output from 'simulate_VAR'.
#' @param p integer; time series lag.
#' @export
VAR_model_matrix <- function(data0, p) {
  eff_T <- ncol(data0) - p
  n <- nrow(data0)
  I_n <- diag(n)
  purrr::map2(
    .x = seq(eff_T), .y = p + seq(eff_T) - 1,
    .f = ~ I_n %x% t( c(1, matrixcalc::vec(data0[, .y:.x, drop = F])) )
  ) %>%
    do.call(rbind, .)
}


# Inference
#' Model inference for Bayesian Vector-AutoRegressive (VAR) model with
#' normal-inverse-Wishart priors.
#' prior for the mean and inverse gamma for the variance.
#' @param data0 A numeric matrix; the multivariate time series. The number
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
VAR_Gibbs <- function(data0, p, b_0, V, v_0, S_0, init_Sigma,
                      num_steps = 1e4, burn_ins = 1e3) {
  if (missing(init_Sigma))
    init_Sigma <- MCMCpack::riwish(v_0, S_0)

  # Initialisation
  T0 <- ncol(data0) - p
  n <- nrow(data0)
  y <- matrixcalc::vec(data0[,seq(T0) + p])
  X <- VAR_model_matrix(data0, p)
  inv_V <- solve(V)
  I_T <- diag(T0)
  Sigma_g <- init_Sigma
  inv_V_times_b_0 <- inv_V %*% b_0
  v_1 <- v_0 + T0

  keep <- num_steps - burn_ins
  res <- vector("list", num_steps)

  pb <- txtProgressBar(1, num_steps, style = 3)
  for (i in 1:num_steps) {
    # Update beta
    inv_Sigma_g <- solve(Sigma_g)
    tX_fac <- t(X) %*% (I_T %x% inv_Sigma_g)
    K_g <- inv_V + tX_fac %*% X
    inv_K_g <- solve(K_g)
    b_g <- inv_K_g %*% (inv_V_times_b_0 + tX_fac %*% y)
    beta_g <- MASS::mvrnorm(1, b_g, inv_K_g)

    # Update Sigma
    coeff_g <- VAR_vec_to_model_coeff(beta_g, n)
    print(microbenchmark::microbenchmark(
      residuals_cov2(y, X, beta_g, n),
      residuals_cov(data0, coeff_g$b_0, coeff_g$B)
    ))
    print(residuals_cov2(y, X, beta_g, n))
    print(residuals_cov(data0, coeff_g$b_0, coeff_g$B))
    co <- readline("Continue?")
    if (co == 'n') return(NULL)
    S_g <- S_0 + residuals_cov(data0, coeff_g$b_0, coeff_g$B)
    Sigma_g <- MCMCpack::riwish(v_1, S_g)

    # Keep track
    res[[i]] <- list(beta = beta_g, Sigma = Sigma_g)
    setTxtProgressBar(pb, i)
  }

  # Tidy format
  tidy_format <- . %>% tail(keep) %>% do.call(rbind, .)
  list(
    Sigma = res %>% purrr::map(~as.numeric(.x$Sigma)) %>% tidy_format(),
    beta = res %>% purrr::map(~.x$beta) %>% tidy_format()
  )
}


residuals_cov <- function(y, X, beta_g, n) {
  tmp <- matrix(y - X %*% beta_g, nrow = n, byrow = F)
  res <- 0
  for (i in 1:ncol(tmp)) {
    res <- res + tmp[,i] %*% t(tmp[,i])
  }
  res
}


#' Compute the posterior mean of the parameters
#' @param model0 Output from 'VAR_Gibbs'.
#' @param n integer; dimension of the data.
#' @param p integer; lag of the time series model.
#' @export
param_posterior_mean <- function(model0, n, p) {
  res <- VAR_vec_to_model_coeff(apply(model0$beta, 2 , mean), p)
  res$Sigma <- matrix(apply(model0$Sigma, 2 , mean), n, n)
  res
}
