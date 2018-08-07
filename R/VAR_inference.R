#' #' Model inference for Bayesian Vector-AutoRegressive (VAR) model with
#' #' independent normal-inverse-Wishart priors.
#' #' @param data0 Matrix; each row is one observation, each column is one measurement / predictor.
#' #' @param lag Integer; the lag of the time series model.
#' #' @param b_0 A numeric vector; the mean for the multivariate normal prior.
#' #' @param B_0 A numeric matrix; the covariance for the multivariate normal prior.
#' #' @param v_0 scalar; degree of freedom for the inverse-Wishart prior.
#' #' @param S_0 matrix; the scale matrix for the inverse-Wishart prior.
#' #' @param init_Sigma (Optional) matrix; the starting value of the noise inverse-covariance.
#' #' @param num_steps integer; number of MCMC steps.
#' #' @return A list of two matrices. The first matrix contains the posterior samples
#' #' of beta and the second matrix contains the posterior samples of Sigma. In both
#' #' matrices, each row corresponding to 1 sample of the vectorised parameter.
#' #' @examples
#' #' sim <- simulate_VAR(m = 1000, 3, 2) # 1000 pts of 3-dimensional lag-2 time series
#' #' data0 <- sim$data
#' #'
#' #' dim_data <- 3
#' #' lag <- 2
#' #' d <- dim_data + dim_data ^2 * lag  # dimension of the vectorised parameter
#' #'
#' #' res <- VAR_Gibbs(
#' #'   data0, lag,
#' #'   b_0 = rnorm(d), B_0 = pdmatrix(d)$Sigma,
#' #'   v_0 = 3, S_0 = pdmatrix(dim_data)$Sigma,
#' #'   init_Sigma = diag(dim_data)
#' #' )
#' #' @export
#' VAR_Gibbs <- function(data0, lag, b_0, B_0, v_0, S_0, init_Sigma, num_steps = 3e3) {
#'   data0 <- train_data(data0, lag)
#'   Y <- data0$Y
#'   X <- data0$X
#'
#'   T0 <- nrow(Y)
#'   n <- ncol(Y)
#'   p <- lag
#'
#'   K <- 1 + n * p
#'   d <- n * K
#'
#'   invB <- solve(B_0)
#'   invVb <- invB %*% b_0
#'   v_1 <- v_0 + T0
#'   if (missing(init_Sigma))
#'     init_Sigma <- MCMCpack::rwish(v_0, S_0)
#'   Sigma <- init_Sigma
#'
#'   XTX <- crossprod(X)
#'   XTY <- crossprod(X, Y)
#'   YTY <- crossprod(Y)
#'
#'   total <- num_steps
#'   Z <- matrix(rnorm(d * total), d, total)
#'   Z2 <- array(rnorm(n * n * total), dim = c(n, n, total))
#'   C <- matrix(0, n, total)
#'   for (i in 1:n) {
#'     C[i,] <- rchisq(total, v_1 - i + 1)
#'   }
#'   C <- sqrt(C)
#'   res_beta <- matrix(0, total, length(b_0))
#'   res_sigma <- matrix(0, total, length(init_Sigma))
#'
#'   get_lower_tri <- . %>% {.[upper.tri(., T)] <- 0; .}
#'
#'   pb <- txtProgressBar(1, num_steps, style = 3)
#'   for (i in 1:total) {
#'     # Update beta
#'     Vg <- solve(invB + kronecker(Sigma, XTX))
#'     bg <- Vg %*% (invVb + as.numeric(XTY %*% Sigma))
#'     beta_g <- bg + t(chol(Vg)) %*% Z[,i]
#'
#'     # Update Sigma
#'     b <- matrix(beta_g, K, n)
#'     S <- S_0 + YTY - t(b) %*% XTY- t(t(b) %*% XTY) + t(b) %*% XTX %*% b
#'     L <- t(chol(solve(S)))
#'     A <- diag(C[, i]) + get_lower_tri(Z2[,,i])
#'     LA <- L %*% A
#'     Sigma <- LA %*% t(LA)
#'
#'     # Keep track
#'     res_beta[i, ] <- beta_g
#'     res_sigma[i, ] <- as.numeric(solve(Sigma))
#'     setTxtProgressBar(pb, i)
#'   }
#'   # Return
#'   list(beta = res_beta, Sigma = res_sigma)
#' }
