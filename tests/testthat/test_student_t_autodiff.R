student_t_test_fun <- function(X, y, b_0, B_0, alpha_0, delta_0, nu,
                               init_beta, init_sigma, num_steps = 1e4) {
  if (missing(init_sigma))
    init_sigma <- 1 / sqrt(rgamma(1, alpha_0 / 2, delta_0 / 2))
  if (missing(init_beta))
    init_beta <- MASS::mvrnorm(1, b_0, B_0)

  # Initialisation
  n <- length(y)
  alpha_1 <- alpha_0 + n
  nu_1 <- nu + 1
  beta_g <- init_beta
  sigma_g <- init_sigma
  inv_B_0 <- solve(B_0)
  inv_B_0_times_b_0 <- inv_B_0 %*% b_0
  beta_res <- matrix(0, num_steps, length(b_0))
  sigma_res <- numeric(num_steps)

  pb <- txtProgressBar(1, num_steps, style = 3)
  for (i in 1:num_steps) {
    #Update lambda
    nu_2 <- nu + sigma_g^(-2) * (y - X %*% beta_g) ^ 2
    Lambda <- rgamma(n, nu_1 / 2, nu_2 / 2)

    # Update beta
    XTL <- A_times_diag_v0(t(X), Lambda)
    XTLX <- XTL %*% X
    XTLy <- XTL %*% y
    B_g <- solve(sigma_g^(-2) * XTLX + inv_B_0)
    b_g <- B_g %*% (sigma_g^(-2) * XTLy + inv_B_0_times_b_0)
    # beta_g <- MASS::mvrnorm(1, b_g, B_g)
    beta_g <- b_g + t(chol(B_g)) %*% rnorm(length(b_g))

    # Update sigma
    delta_g <- delta_0 + A_times_diag_v0(t(y - X %*% beta_g), Lambda) %*% (y - X %*% beta_g)
    sigma_g <- 1 / sqrt(rgamma(1, alpha_1 / 2, delta_g / 2))

    # Keep track
    beta_res[i,] <- beta_g
    sigma_res[i] <- sigma_g
    setTxtProgressBar(pb, i)
  }

  # Tidy format
  list(sigma = sigma_res, beta = beta_res)
}


testthat::test_that("Testing student-t autodiff", {
  skip_on_cran()
  set.seed(123);
  n <- 1000
  p <- 5
  data0 <- student_t_data(n, p, intercept = FALSE)

  set.seed(123);
  b_0 <- rnorm(p)
  B_0 <- pdmatrix(p)$Sigma
  h <- 0.000001

  res <- student_t_AD(
    data0$X, data0$y, b_0 = b_0, B_0 = B_0,
    alpha_0 = 3, delta_0 = 5, nu = 20, num_steps = 200, burn_ins = 0
  )

  library(magrittr)
  tailMeans <- . %>% tail(1800) %>% colMeans()
  tailmean <- . %>% tail(1800) %>% mean()
  #================== Check sensitivity of beta and sigma wrt b0 ======================
  for (i in 1:p) {
    new_b_0 <- b_0
    new_b_0[i] <- new_b_0[i] + h

    set.seed(123);
    a <- student_t_test_fun(
      data0$X, data0$y, b_0 = b_0, B_0 = B_0,
      alpha_0 = 3, delta_0 = 5, nu = 20, num_steps = 2000
    )

    set.seed(123);
    b <- student_t_test_fun(
      data0$X, data0$y, b_0 = new_b_0, B_0 = B_0,
      alpha_0 = 3, delta_0 = 5, nu = 20, num_steps = 2000
    )

    v1 <- (tailMeans(b$beta) - tailMeans(a$beta)) / h
    v2 <- matrix(colMeans(res$d_beta$d_b0), p, p)[,i]
    testthat::expect_lt(sum(abs(v1 - v2)), 1e-5)

    s1 <- (tailmean(b$sigma^2) - tailmean(a$sigma^2)) / h
    s2 <- colMeans(res$d_sigma2$d_b0)[i]
    testthat::expect_lt(sum(abs(s1 - s2)), 1e-5)

    # print(rbind(v1, v2))
    # print(rbind(s1, s2))
  }

  #================== Check sensitivity of beta and sigma wrt B_0 =====================
  for (j in 1:p) {
    for (i in 1:p) {
      new_B_0 <- B_0
      new_B_0[i, j] <- new_B_0[i, j] + h

      set.seed(123);
      a <- student_t_test_fun(
        data0$X, data0$y, b_0 = b_0, B_0 = B_0,
        alpha_0 = 3, delta_0 = 5, nu = 20, num_steps = 2000
      )

      set.seed(123);
      b <- student_t_test_fun(
        data0$X, data0$y, b_0 = b_0, B_0 = new_B_0,
        alpha_0 = 3, delta_0 = 5, nu = 20, num_steps = 2000
      )

      v1 <- (tailMeans(b$beta) - tailMeans(a$beta)) / h
      v2 <- matrix(colMeans(res$d_beta$d_B0), p, p^2)[, i + (j-1) * p]
      testthat::expect_lt(sum(abs(v1 - v2)), 1e-5)

      s1 <- (tailmean(b$sigma^2) - tailmean(a$sigma^2)) / h
      s2 <- colMeans(res$d_sigma2$d_B0)[i + (j-1) * p]
      testthat::expect_lt(sum(abs(s1 - s2)), 1e-5)

      # print(rbind(v1, v2))
      # print(rbind(s1, s2))
    }
  }

  #================== Check sensitivity of beta and sigma wrt alpha_0 =================
  set.seed(123);
  a <- student_t_test_fun(
    data0$X, data0$y, b_0 = b_0, B_0 = B_0,
    alpha_0 = 3, delta_0 = 5, nu = 20, num_steps = 5000
  )

  set.seed(123);
  b <- student_t_test_fun(
    data0$X, data0$y, b_0 = b_0, B_0 = B_0,
    alpha_0 = 3 + h, delta_0 = 5, nu = 20, num_steps = 5000
  )

  v1 <- (colMeans(tail(b$beta, 4e3)) - colMeans(tail(a$beta, 4e3))) / h
  v2 <- colMeans(res$d_beta$d_alpha0)
  testthat::expect_lt(sum(abs(v1 - v2)), 1e-5)

  s1 <- (tailmean(b$sigma^2) - tailmean(a$sigma^2)) / h
  s2 <- colMeans(res$d_sigma2$d_alpha0)
  testthat::expect_lt(sum(abs(s1 - s2)), 1e-5)

  # print(rbind(v1, v2))
  # print(rbind(s1, s2))

  #================== Check sensitivity of beta and sigma wrt delta_0 =================
  set.seed(123);
  a <- student_t_test_fun(
    data0$X, data0$y, b_0 = b_0, B_0 = B_0,
    alpha_0 = 3, delta_0 = 5, nu = 20, num_steps = 5000
  )

  set.seed(123);
  b <- student_t_test_fun(
    data0$X, data0$y, b_0 = b_0, B_0 = new_B_0,
    alpha_0 = 3, delta_0 = 5 + h, nu = 20, num_steps = 5000
  )

  v1 <- (colMeans(tail(b$beta, 4e3)) - colMeans(tail(a$beta, 4e3))) / h
  v2 <- colMeans(res$d_beta$d_delta0)
  testthat::expect_lt(sum(abs(v1 - v2)), 1e-5)

  s1 <- (tailmean(b$sigma^2) - tailmean(a$sigma^2)) / h
  s2 <- colMeans(res$d_sigma2$d_delta0)
  testthat::expect_lt(sum(abs(s1 - s2)), 1e-5)

  # print(rbind(v1, v2))
  # print(rbind(s1, s2))
})
