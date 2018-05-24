testthat::test_that("Testing Wishart autodiff", {
  skip_cran()
  set.seed(123)
  n <- 1000
  p <- 2
  k <- 3
  beta <- rnorm(p + 2)
  gamma <- rnorm(k + 1)
  sigma <- matrix(c(1, 0.2, 0.2, 1), 2, 2)
  data0 <- wishart_data(n, p, k, beta = beta, gamma = gamma, sigma = sigma,
                        intercept_1 = TRUE, intercept_2 = TRUE)

  set.seed(234)
  b_0 <- rnorm(p + 2)
  B_0 <- pdmatrix(p+2)$Sigma
  g_0 <- rnorm(k + 1)
  G_0 <- pdmatrix(k + 1)$Sigma
  v_0 <- 5
  R_0 <- pdmatrix(2)$Sigma
  seed <- 1234

  set.seed(seed)
  res <- wishart_test_fun(Xy = data0$Xy, y = data0$y, Xs = data0$Xs, s = data0$s,
                       b_0 = b_0, B_0 = B_0, g_0 = g_0, G_0 = G_0,
                       v_0 = v_0, R_0 = R_0)

  res_AD <- wishart_AD(Xy = data0$Xy, y = data0$y, Xs = data0$Xs, s = data0$s,
                       b_0 = b_0, B_0 = B_0,
                       g_0 = g_0, G_0 = G_0,
                       v_0 = v_0, R_0 = R_0,
                       num_steps = 1e3, burn_ins = 0)

  h <- 1e-5
  for (i in 1:length(b_0)) {
    new_b0 <- b_0
    new_b0[i] <- new_b0[i] + h
    set.seed(seed)
    res2 <- wishart_test_fun(Xy = data0$Xy, y = data0$y, Xs = data0$Xs, s = data0$s,
                         b_0 = new_b0, B_0 = B_0, g_0 = g_0, G_0 = G_0,
                         v_0 = v_0, R_0 = R_0)

    fd_sen <- (colMeans(res2$beta) - colMeans(res$beta)) /  h
    ad_sen <- res_AD$d_beta$d_b0

    s1 <- (colMeans(res2$sigma) - colMeans(res$sigma)) / h
    s2 <- matrix(colMeans(res_AD$d_Sigma$d_b0), 4, 4)[,i]

    print(rbind(fd_sen, ad_sen = matrix(colMeans(ad_sen), 4, 4)[,i]))
    print(rbind(s1, s2))
  }
})
