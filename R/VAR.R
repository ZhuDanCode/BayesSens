#' Simulate data from a Vector AutoRegressive (VAR) model
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


# Inference
#' Model inference for Bayesian Vector-AutoRegressive (VAR) model with
#' normal-inverse-Wishart priors.
#' @param Y Matrix; each row is one observation, each column is one measurement / predictor.
#' @param X Matrix;
#' @param lag Integer; the lag of the time series model.
#' @param b_0 Vector; the intercept
#' @param b_0 A numeric vector; the mean for the multivariate normal prior.
#' @param B_0 A numeric matrix; the covariance for the multivariate normal prior.
#' @param v_0 scalar; degree of freedom for the inverse-Wishart prior.
#' @param S_0 matrix; the scale matrix for the inverse-Wishart prior.
#' @param init_Sigma (Optional) matrix; the starting value of the noise covariance.
#' @param num_steps integer; number of MCMC steps.
#' @param burn_ins integer; number of burn-ins.
VAR_Gibbs_2 <- function(Y, X, lag, b_0, B_0, v_0, S_0, init_Sigma,
                        num_steps = 3e3, burn_ins = 1e3) {
  T0 <- nrow(Y)
  n <- ncol(Y)
  p <- lag

  K <- 1 + n * p
  d <- n * K

  invB <- solve(B_0)
  invVb <- invB %*% b_0
  v_1 <- v_0 + T0
  Sigma <- init_Sigma

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

  pb <- txtProgressBar(1, num_steps, style = 3)
  for (i in 1:total) {
    # Update beta
    Vg <- solve(invB + kronecker(Sigma, XTX))
    bg <- Vg %*% (invVb + as.numeric(XTY %*% Sigma))
    beta_g <- bg + t(chol(Vg)) %*% Z[,i]

    # Update Sigma
    b <- matrix(beta_g, K, n)
    S <- S_0 + YTY - t(b) %*% XTY- t(t(b) %*% XTY) + t(b) %*% XTX %*% b
    L <- t(chol(solve(S)))
    A <- diag(C[, i]) + get_lower_tri(Z2[,,i])
    LA <- L %*% A
    Sigma <- LA %*% t(LA)

    # Keep track
    res_beta[i, ] <- beta_g
    res_sigma[i, ] <- as.numeric(solve(Sigma))
    setTxtProgressBar(pb, i)
  }
  # Return
  list(beta = res_beta, Sigma = res_sigma)
}

train_data <- function(Y, lag) {
  nr <- nrow(Y)
  nc <- ncol(Y)
  X <- matrix(0, nr - lag, lag * nc)
  for (i in 1:lag) {
    X[, (1 + (i-1) * nc):(i * nc)] <- Y[(1 + lag - i):(nr - i), ]
  }
  list(X = cbind(1, X), Y = Y[(1 + lag):nr, ])
}

# Example / Unit test
# X <- read.csv("../Bayesian MCMC Sensitivity/var paper/code/X.csv", header = F)
# Y <- read.csv("../Bayesian MCMC Sensitivity/var paper/code/Y.csv", header = F)
# Y <- matrix(unlist(Y), ncol = 3, byrow = T)
# X <- cbind(1, as.matrix(X))
#
# n <- ncol(Y)
# p <- 2
# adj_dim <- n * (1 + n * p)
# b_0 <- rnorm(adj_dim)
# B_0 <- 0.1 * diag(adj_dim)
# v_0 <- 3
# S_0 <- 0.1 * diag(n)
# init_Sigma <- 0.1 * diag(n)
#
# res <- VAR_Gibbs_2(Y, X, lag = 2, b_0, B_0, v_0, S_0, init_Sigma)
#
# sim <- simulate_VAR(m = 1000, n = 3, p = 2)
# data0 <- train_data(t(sim$data), 2)
# res <- VAR_Gibbs_2(data0$Y, data0$X, lag = 2, b_0, B_0, v_0, S_0, init_Sigma)
#
