# This function is needed to align with the implementation in the Autodiff function.
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
  skip_on_cran()
  set.seed(123);
  n <- 1000
  p <- 5
  data0 <- gaussian_data(n, p, intercept = FALSE)

  set.seed(123);
  b_0 <- rnorm(p)
  B_0 <- pdmatrix(p)$Sigma
  h <- 0.00001

  res <- gaussian_AD(
    data0$X, data0$y, b_0 = b_0, B_0 = B_0,
    alpha_0 = 3, delta_0 = 5, num_steps = 200, burn_ins = 0
  )

  #================== Check sensitivity of beta and sigma wrt b0 ======================
  for (i in 1:p) {
    new_b_0 <- b_0
    new_b_0[i] <- new_b_0[i] + h

    set.seed(123);
    a <- test_fun(
      data0$X, data0$y, b_0 = b_0, B_0 = B_0,
      alpha_0 = 3, delta_0 = 5, num_steps = 2000
    )

    set.seed(123);
    b <- test_fun(
      data0$X, data0$y, b_0 = new_b_0, B_0 = B_0,
      alpha_0 = 3, delta_0 = 5, num_steps = 2000
    )

    v1 <- (colMeans(b$beta) - colMeans(a$beta)) / h
    v2 <- matrix(colMeans(res$d_beta$d_b0), p, p)[,i]
    testthat::expect_lt(sum(abs(v1 - v2)), 1e-5)

    s1 <- (mean(b$sigma^2) - mean(a$sigma^2)) / h
    s2 <- colMeans(res$d_sigma2$d_b0)[i]
    testthat::expect_lt(sum(abs(s1 - s2)), 1e-5)
  }

  #================== Check sensitivity of beta and sigma wrt B_0 =====================
  for (j in 1:p) {
    for (i in 1:p) {
      new_B_0 <- B_0
      new_B_0[i, j] <- new_B_0[i, j] + h

      set.seed(123);
      a <- test_fun(
        data0$X, data0$y, b_0 = b_0, B_0 = B_0,
        alpha_0 = 3, delta_0 = 5, num_steps = 2000
      )

      set.seed(123);
      b <- test_fun(
        data0$X, data0$y, b_0 = b_0, B_0 = new_B_0,
        alpha_0 = 3, delta_0 = 5, num_steps = 2000
      )

      v1 <- (colMeans(b$beta) - colMeans(a$beta)) / h
      v2 <- matrix(colMeans(res$d_beta$d_B0), p, p^2)[, i + (j-1) * p]
      testthat::expect_lt(sum(abs(v1 - v2)), 1e-5)

      s1 <- (mean(b$sigma^2) - mean(a$sigma^2)) / h
      s2 <- colMeans(res$d_sigma2$d_B0)[i + (j-1) * p]
      testthat::expect_lt(sum(abs(s1 - s2)), 1e-5)
    }
  }

  #================== Check sensitivity of beta and sigma wrt alpha_0 =================
  set.seed(123);
  a <- test_fun(
    data0$X, data0$y, b_0 = b_0, B_0 = B_0,
    alpha_0 = 3, delta_0 = 5, num_steps = 5000
  )

  set.seed(123);
  b <- test_fun(
    data0$X, data0$y, b_0 = b_0, B_0 = B_0,
    alpha_0 = 3 + h, delta_0 = 5, num_steps = 5000
  )

  v1 <- (colMeans(tail(b$beta, 4e3)) - colMeans(tail(a$beta, 4e3))) / h
  v2 <- colMeans(res$d_beta$d_alpha0)
  testthat::expect_lt(sum(abs(v1 - v2)), 1e-5)

  s1 <- (mean(b$sigma^2) - mean(a$sigma^2)) / h
  s2 <- colMeans(res$d_sigma2$d_alpha0)
  testthat::expect_lt(sum(abs(s1 - s2)), 1e-5)

  #================== Check sensitivity of beta and sigma wrt delta_0 =================
  set.seed(123);
  a <- test_fun(
    data0$X, data0$y, b_0 = b_0, B_0 = B_0,
    alpha_0 = 3, delta_0 = 5, num_steps = 5000
  )

  set.seed(123);
  b <- test_fun(
    data0$X, data0$y, b_0 = b_0, B_0 = new_B_0,
    alpha_0 = 3, delta_0 = 5 + h, num_steps = 5000
  )

  v1 <- (colMeans(tail(b$beta, 4e3)) - colMeans(tail(a$beta, 4e3))) / h
  v2 <- colMeans(res$d_beta$d_delta0)
  testthat::expect_lt(sum(abs(v1 - v2)), 1e-5)

  s1 <- (mean(b$sigma^2) - mean(a$sigma^2)) / h
  s2 <- colMeans(res$d_sigma2$d_delta0)
  testthat::expect_lt(sum(abs(s1 - s2)), 1e-5)
})
