context("wishart: Check Auto-differentiation results are consistent with Numerical differentiation")
library(BayesSense)

# Wishart Gibbs with Wishart draws implemented via Bartlett's decomposition.
# This function is needed to align with the implementation in the Autodiff function.
test_fun <- function(Xy, y, b_0, B_0, Xs, s, g_0, G_0, v_0, R_0,
                     init_gamma, init_Sigma, num_steps = 1e4) {
  if (missing(init_gamma))
    init_gamma <- g_0 + t(chol(G_0)) %*% rnorm(length(g_0)) %>% as.vector()
  if (missing(init_Sigma))
    init_Sigma <- matrix(rWishart(1, v_0, R_0), 2, 2)

  # Initialisation
  n <- length(y)
  len_b0 <- length(b_0)
  len_g0 <- length(g_0)
  v_1 <- v_0 + n
  gamma_g <- init_gamma
  sigma_g <- init_Sigma
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
    w_11_g <- sigma_g[1] - sigma_g[2] * sigma_g[3] / sigma_g[4]
    b_sum <- XyTy - sigma_g[2] * sigma_g[3] / sigma_g[4] * (XyTs - XyTXs %*% gamma_g)
    B_g <- solve(XyTXy / w_11_g + inv_B_0)
    b_g <- B_g %*% (b_sum / w_11_g + inv_B_0_times_b_0)
    chol_Bg <- t(chol(B_g))
    beta_g <- b_g + chol_Bg %*% Zy[, i]

    # Update gamma
    w_22_g <- sigma_g[4] - sigma_g[2] * sigma_g[3] / sigma_g[1]
    g_sum <- XsTs - sigma_g[2] * sigma_g[3] / sigma_g[1] * (XsTy - XsTXy %*% beta_g)
    G_g <- solve(XsTXs / w_22_g + inv_G_0)
    g_g <- G_g %*% (g_sum / w_22_g + inv_G_0_times_g_0)
    chol_Gg <- t(chol(G_g))
    gamma_g <- g_g + chol_Gg %*% Zs[, i]

    # Update sigma
    delta <- matrix(c(beta_g, rep(0, len_b0 + len_g0), gamma_g), ncol = 2)
    vTX_delta <- vTX %*% delta
    X_delta <- X %*% delta
    R_1 <- solve(inv_R_0 + vTv - vTX_delta - t(vTX_delta) + crossprod(X_delta))
    L <- t(chol(R_1))
    A <- diag(C[, i]) + get_lower_tri(Z2[,,i])
    LA <- L %*% A
    sigma_g <- solve(LA %*% t(LA))

    # Keep track
    res[[i]] <- list(beta = beta_g, gamma = gamma_g, Sigma = sigma_g)
    setTxtProgressBar(pb, i)
  }

  # Tidy format
  extract_from <- function(a, b) {
    map_reduce(b, ~.x %>% extract2(a) %>% as.vector(), rbind)
  }
  map_named(c("beta", "gamma", "Sigma"), ~extract_from(.x, res))
}

# num_diff_pipe :: list <configuration> -> list <matrices>
num_diff_pipe <- function(config_list, base_case, h) {
  config_list %>%
    purrr::map(~num_diffs(.x, base_case, h)) %>%
    tidy_list() %>%
    purrr::map(t)
}

# num_diffs :: list -> Obj -> list <summary statistics>
num_diffs <- function(new_parameters, base_case, h) {
  set.seed(345)
  new_case <- do.call(test_fun, new_parameters)
  list(d_beta = num_diff(new_case$beta, base_case$beta, h, 0.9),
       d_gamma = num_diff(new_case$gamma, base_case$gamma, h, 0.9),
       d_Sigma = num_diff(new_case$Sigma, base_case$Sigma, h, 0.9))
}

# Test if numerical differentiation and auto-differentiation give
# a difference less than a threshold.
test_diff <- function(num_diff, auto_diff) {
  x <- unlist(num_diff)
  y <- unlist(auto_diff)
  print(cbind(num_diff = x, auto_diff = y))
  print(max(abs(x - y)))
  expect_true(all(abs(x - y) < 1e-5))
}

testthat::test_that("Testing Wishart autodiff", {
  # skip_on_cran()
  skip_if_not(Sys.getenv("run_all_tests", TRUE))
  #======================== Generate data ==========================================
  set.seed(345)
  n <- 1000
  p <- 2
  k <- 3
  beta <- rnorm(p + 2)
  gamma <- rnorm(k + 1)
  Sigma <- matrix(c(1, 0.2, 0.2, 1), 2, 2)
  data0 <- wishart_data(n, p, k, beta = beta, gamma = gamma, Sigma = Sigma,
                        intercept_1 = TRUE, intercept_2 = TRUE)

  #======================= Generate parameters =====================================
  set.seed(234)
  b_0 <- numeric(p + 2)
  B_0 <- diag(p + 2)
  g_0 <- numeric(k + 1)
  G_0 <- diag(k + 1)
  v_0 <- 5
  R_0 <- diag(2)
  h <- 1e-7
  default_parameters <- list(Xy = data0$Xy, y = data0$y, Xs = data0$Xs, s = data0$s,
                        b_0 = b_0, B_0 = B_0, g_0 = g_0, G_0 = G_0,
                        v_0 = v_0, R_0 = R_0)
  seed <- 345
  set.seed(seed)
  base_case <- do.call(test_fun, default_parameters)

  set.seed(seed)
  res <- do.call(wishart_AD, append(default_parameters, list(num_steps = 2e3)))

  #================ Check sensitivity of beta, gamma and Sigma wrt b0 ==================
  print("Sensitivity of beta, gamma, Sigma wrt b0")
  num_diff_sense <- 1:length(b_0) %>%
    # Generate Configurations
    purrr::map(function(i) {
      new_parameters <- default_parameters
      new_parameters$b_0[i] <- new_parameters$b_0[i] + h
      new_parameters
    }) %>%
    num_diff_pipe(base_case, h)
  auto_diff_sense <- list(
    d_beta = auto_diff(res, "d_beta", "d_b0"),
    d_gamma = auto_diff(res, "d_gamma", "d_b0"),
    d_Sigma = auto_diff(res, "d_Sigma", "d_b0")
  )
  test_diff(num_diff_sense, auto_diff_sense)

  #=============== Check sensitivity of beta, gamma and Sigma wrt B0 ======================
  print("Sensitivity of beta, gamma, Sigma wrt B0")
  num_diff_sense <- expand.grid(1:nrow(B_0), 1:ncol(B_0)) %>%
    {purrr::map2(.[,1], .[,2], function(i, j) {
      new_parameters <- default_parameters
      new_parameters$B_0[i, j] <- new_parameters$B_0[i, j] + h
      if (i != j) {
        new_parameters$B_0[j, i] <- new_parameters$B_0[j, i] + h
      }
      new_parameters
    })} %>%
    num_diff_pipe(base_case, h)
  auto_diff_sense <- list(
    d_beta = auto_diff(res, "d_beta", "d_B0"),
    d_gamma = auto_diff(res, "d_gamma", "d_B0"),
    d_Sigma = auto_diff(res, "d_Sigma", "d_B0")
  )
  test_diff(num_diff_sense, auto_diff_sense)

  #=============== Check sensitivity of beta, gamma and Sigma wrt g0 ======================
  print("Sensitivity of beta, gamma, Sigma wrt g0")
  num_diff_sense <- 1:(k+1) %>%
    # Generate Configurations
    purrr::map(function(i) {
      new_parameters <- default_parameters
      new_parameters$g_0[i] <- new_parameters$g_0[i] + h
      new_parameters
    }) %>%
    num_diff_pipe(base_case, h)
  auto_diff_sense <- list(
    d_beta = auto_diff(res, "d_beta", "d_g0"),
    d_gamma = auto_diff(res, "d_gamma", "d_g0"),
    d_Sigma = auto_diff(res, "d_Sigma", "d_g0")
  )
  test_diff(num_diff_sense, auto_diff_sense)

  #=============== Check sensitivity of beta, gamma and Sigma wrt G0 ======================
  print("Sensitivity of beta, gamma, Sigma wrt G_0")
  num_diff_sense <- expand.grid(1:(k+1), 1:(k+1)) %>%
    {purrr::map2(.[,1], .[,2], function(i, j) {
      new_parameters <- default_parameters
      new_parameters$G_0[i, j] <- new_parameters$G_0[i, j] + h
      if (i != j) {
        new_parameters$G_0[j, i] <- new_parameters$G_0[j, i] + h
      }
      new_parameters
    })} %>%
    num_diff_pipe(base_case, h)
  auto_diff_sense <- list(
    d_beta = auto_diff(res, "d_beta", "d_G0"),
    d_gamma = auto_diff(res, "d_gamma", "d_G0"),
    d_Sigma = auto_diff(res, "d_Sigma", "d_G0")
  )
  test_diff(num_diff_sense, auto_diff_sense)

  #=============== Check sensitivity of beta, gamma and Sigma wrt v0 ======================
  print("Sensitivity of beta, gamma, Sigma wrt v_0")
  new_parameters <- default_parameters
  new_parameters$v_0 <- new_parameters$v_0 + h

  num_diff_sense <- list(new_parameters) %>% num_diff_pipe(base_case, h)

  auto_diff_sense <- list(
    d_beta = auto_diff(res, "d_beta", "d_v0"),
    d_gamma = auto_diff(res, "d_gamma", "d_v0"),
    d_Sigma = auto_diff(res, "d_Sigma", "d_v0")
  )
  test_diff(num_diff_sense, auto_diff_sense)

  #=============== Check sensitivity of beta, gamma and Sigma wrt R0 ======================
  print("Sensitivity of beta, gamma, Sigma wrt R0")
  num_diff_sense <- expand.grid(1:2, 1:2) %>%
    {purrr::map2(.[,1], .[,2], function(i, j) {
      new_parameters <- default_parameters
      new_parameters$R_0[i, j] <- new_parameters$R_0[i, j] + h
      if (i != j) {
        new_parameters$R_0[j, i] <- new_parameters$R_0[j, i] + h
      }
      new_parameters
    })} %>%
    num_diff_pipe(base_case, h)
  auto_diff_sense <- list(
    d_beta = auto_diff(res, "d_beta", "d_R0"),
    d_gamma = auto_diff(res, "d_gamma", "d_R0"),
    d_Sigma = auto_diff(res, "d_Sigma", "d_R0")
  )
  test_diff(num_diff_sense, auto_diff_sense)

})
