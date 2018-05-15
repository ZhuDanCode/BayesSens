# Implementation two
testthat::test_that("Test gaussian autodiff", {
  set.seed(123)
  n <- 1000
  p <- 5
  data0 <- gaussian_data(n, p, intercept = TRUE)

  set.seed(123)
  res <- gaussian_AD(data0$X, data0$y,
    b_0 = rnorm(p+1), B_0 = pdmatrix(p+1)$Sigma,  # add one for the intercept
    alpha_0 = 13, delta_0 = 8, num_steps = 200, burn_ins = 0
  )
  set.seed(123)
  res_2 <- BayesSense:::gaussian_AD_2(data0$X, data0$y,
    b_0 = rnorm(p+1), B_0 = pdmatrix(p+1)$Sigma,  # add one for the intercept
    alpha_0 = 13, delta_0 = 8, num_steps = 200, burn_ins = 0
  )

  for (i in names(res$d_beta)) {
    v1 <- colMeans(res$d_beta[[i]])
    v2 <- colMeans(res_2$d_beta[[i]])
    expect_lt(sum(abs(v1 - v2)), 1e-10)
  }

  for (i in names(res$d_sigma2)) {
    v1 <- colMeans(res$d_sigma2[[i]])
    v2 <- colMeans(res_2$d_sigma2[[i]])
    expect_lt(sum(abs(v1 - v2)), 1e-10)
  }
})
