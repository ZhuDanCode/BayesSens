context("Check Auto-differentiation results are consistent with Numerical differentiation")
library(BayesSense)

# This function aligns with the implementation in the Autodiff function, and it is needed
# for checking against the numerical differentiation.
test_fun <- function(X, y, b_0, B_0, alpha_0, delta_0, init_sigma, num_steps) {
  if (missing(init_sigma))
    init_sigma <- 1 / sqrt(rgamma(1, alpha_0 / 2, delta_0 / 2))

  # Initialisation
  n <- length(y)
  alpha_1 <- alpha_0 + n
  sigma_g <- init_sigma
  inv_B_0 <- solve(B_0)
  inv_B_0_times_b_0 <- inv_B_0 %*% b_0
  XTX <- crossprod(X)
  XTy <- crossprod(X, y)
  beta_res <- matrix(0, num_steps, length(b_0))
  sigma_res <- numeric(num_steps)

  pb <- txtProgressBar(1, num_steps, style = 3)
  for (i in 1:num_steps) {
    # Update beta
    B_g <- solve(sigma_g^(-2) * XTX + inv_B_0)
    b_g <- B_g %*% (sigma_g^(-2) * XTy + inv_B_0_times_b_0)
    # beta_g <- MASS::mvrnorm(1, b_g, B_g)
    z <- rnorm(length(b_0))
    beta_g <- b_g + t(chol(B_g)) %*% z

    # Update sigma
    delta_g <- delta_0 + sum((y - X %*% beta_g)^2)
    # sigma_g <- 1 / sqrt(rgamma(1, alpha_1 / 2, delta_g / 2))
    G <- rgamma(1, alpha_1 / 2, 1)
    sigma_g <- 1 / sqrt(2 / delta_g * G)

    # Keep track
    beta_res[i,] <- beta_g
    sigma_res[i] <- sigma_g
    setTxtProgressBar(pb, i)
  }

  list(sigma = sigma_res, beta = beta_res)
}

testthat::test_that("Testing gaussian autodiff", {
  # testthat::skip_on_cran()
  #======================= Setup and Initialisation  ==============================
  library(BayesSense)
  # Setup data
  set.seed(123)
  n <- 1000
  p <- 5
  data0 <- gaussian_data(n, p, intercept = FALSE)

  # Initialise Base Case and Autodiff Case
  set.seed(123)
  b_0 <- rnorm(p)
  B_0 <- pdmatrix(p)$Sigma
  h <- 0.00001
  default_parameters <- list(
    X = data0$X, y = data0$y,
    b_0 = b_0, B_0 = B_0, alpha_0 = 3, delta_0 = 5,
    num_steps = 2000
  )
  set.seed(123)
  base_case <- do.call(test_fun, default_parameters)
  set.seed(123)
  res <- do.call(gaussian_AD, default_parameters)

  # Helper functions
  # This function computes the sensitivity of the posterior mean
  # x , y are matrices with 'num_steps' rows and 'dim(parameters)' cols.
  numdiff <- function(x, y) {
    if (is.matrix(x) && is.matrix(y)) {
      return((colMeans_tail(x) - colMeans_tail(y)) / h)
    }
    (mean_tail(x) - mean_tail(y)) / h
  }

  #================== Check sensitivity of beta and sigma wrt b0 ======================
  for (i in 1:p) {
    new_parameters <- default_parameters
    new_parameters$b_0[i] <- new_parameters$b_0[i] + h

    set.seed(123)
    new_case <- do.call(test_fun, new_parameters)
    v1 <- numdiff(new_case$beta, base_case$beta)
    s1 <- numdiff(new_case$sigma^2, base_case$sigma^2)

    # Compare with autodiff
    v2 <- matrix(colMeans_tail(res$d_beta$d_b0), p, p)[,i]
    testthat::expect_lt(sum(abs(v1 - v2)), 1e-5)
    s2 <- colMeans_tail(res$d_sigma2$d_b0)[i]
    testthat::expect_lt(sum(abs(s1 - s2)), 1e-5)
    print(rbind(numdiff = v1, autodiff = v2))
    print(rbind(numdiff = s1, autodiff = s2))
  }

  #================== Check sensitivity of beta and sigma wrt B_0 =====================
  for (j in 1:p) {
    for (i in 1:p) {
      new_parameters <- default_parameters
      new_parameters$B_0[i, j] <- new_parameters$B_0[i, j] + h

      set.seed(123)
      new_case <- do.call(test_fun, new_parameters)
      v1 <- numdiff(new_case$beta, base_case$beta)
      s1 <- numdiff(new_case$sigma^2, base_case$sigma^2)

      # Compare with autodiff
      v2 <- matrix(colMeans_tail(res$d_beta$d_B0), p, p^2)[, i + (j-1) * p]
      testthat::expect_lt(sum(abs(v1 - v2)), 1e-5)
      s2 <- colMeans_tail(res$d_sigma2$d_B0)[i + (j-1) * p]
      testthat::expect_lt(sum(abs(s1 - s2)), 1e-5)
      print(rbind(numdiff = v1, autodiff = v2))
      print(rbind(numdiff = s1, autodiff = s2))
    }
  }

  #================== Check sensitivity of beta and sigma wrt alpha_0 =================
  default_parameters$num_steps <- 5000
  new_parameters <- default_parameters
  new_parameters$alpha_0 <- new_parameters$alpha_0 + h

  set.seed(123)
  base_case <- do.call(test_fun, default_parameters)
  set.seed(123)
  new_case <- do.call(test_fun, new_parameters)
  v1 <- numdiff(new_case$beta, base_case$beta)
  s1 <- numdiff(new_case$sigma^2, base_case$sigma^2)

  # Compare with autodiff
  v2 <- colMeans_tail(res$d_beta$d_alpha0)
  testthat::expect_lt(sum(abs(v1 - v2)), 1e-5)
  s2 <- colMeans_tail(res$d_sigma2$d_alpha0)
  testthat::expect_lt(sum(abs(s1 - s2)), 1e-5)
  print(rbind(numdiff = v1, autodiff = v2))
  print(rbind(numdiff = s1, autodiff = s2))

  #================== Check sensitivity of beta and sigma wrt delta_0 =================
  default_parameters$num_steps <- 5000
  new_parameters <- default_parameters
  new_parameters$delta_0 <- new_parameters$delta_0 + h

  set.seed(123)
  base_case <- do.call(test_fun, default_parameters)
  set.seed(123)
  new_case <- do.call(test_fun, new_parameters)
  v1 <- numdiff(new_case$beta, base_case$beta)
  s1 <- numdiff(new_case$sigma^2, base_case$sigma^2)

  # Compare with autodiff
  v2 <- colMeans_tail(res$d_beta$d_delta0)
  testthat::expect_lt(sum(abs(v1 - v2)), 1e-5)
  s2 <- colMeans_tail(res$d_sigma2$d_delta0)
  testthat::expect_lt(sum(abs(s1 - s2)), 1e-5)
  print(rbind(numdiff = v1, autodiff = v2))
  print(rbind(numdiff = s1, autodiff = s2))
})
