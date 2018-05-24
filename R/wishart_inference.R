#' Model inference for Wishart regression model with multivariate normal
#' priors for the means and inverse gamma for the variances.
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
#' @param init_sigma (Optional) A 2x2 symmetric matrix; the starting value of the sigma parameter.
#' @param num_steps Integer; number of MCMC steps.
#' @param burn_ins Integer; number of burn-ins.
#' @examples
#' \dontrun{
#' n <- 1000
#' p <- 2
#' k <- 3
#' data0 <- wishart_data(n, p, k,
#'                       sigma = matrix(c(1, 0.2, 0.2, 1), 2, 2),
#'                       intercept_1 = TRUE, intercept_2 = TRUE)
#' res <- wishart_Gibbs(Xy = data0$Xy, y = data0$y,
#'                      b_0 = rnorm(p + 2), B_0 = pdmatrix(p + 2)$Sigma,
#'                      Xs <- data0$Xs, s = data0$s,
#'                      g_0 = rnorm(k + 1), G_0 = pdmatrix(k + 1)$Sigma,
#'                      v_0 = 5, R_0 = pdmatrix(2)$Sigma)
#' }
#' @export
wishart_Gibbs <- function(Xy, y, b_0, B_0, Xs, s, g_0, G_0, v_0, R_0,
                          init_gamma, init_sigma,
                          num_steps = 1e4, burn_ins = 1e3) {
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

  keep <- num_steps - burn_ins
  res <- vector("list", num_steps)

  #pre-simluation
  Zy <- matrix(rnorm(len_b0 * num_steps), len_b0, num_steps)
  Zs <- matrix(rnorm(len_g0 * num_steps), len_g0, num_steps)

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
    R_1 <- solve(inv_R_0 + vTv - vTX_delta - t(vTX_delta) + crossprod(X_delta))
    sigma_g <- solve(rWishart(1, v_1, R_1)[,,1])

    # Keep track
    res[[i]] <- list(beta = beta_g, gamma = gamma_g, sigma = sigma_g)
    setTxtProgressBar(pb, i)
  }

  # Tidy format
  list(
    sigma = res %>% purrr::map(~.x$sigma) %>% tail(keep) %>% lapply(as.vector) %>% do.call(rbind, .),
    beta = res %>% purrr::map(~.x$beta) %>% tail(keep) %>% lapply(as.vector) %>% do.call(rbind, .),
    gamma = res %>% purrr::map(~.x$gamma) %>% tail(keep) %>% lapply(as.vector) %>% do.call(rbind, .)
  )
}
