# Cross-checking the gaussian autodiff using a different implementation
testthat::test_that("Test gaussian autodiff", {
  set.seed(123)
  n <- 1000
  p <- 5
  data0 <- gaussian_data(n, p, intercept = TRUE)

  set.seed(123)
  default_parameters <- list(
    X = data0$X, y = data0$y,
    b_0 = rnorm(p+1), B_0 = pdmatrix(p+1)$Sigma,  # add one for the intercept
    alpha_0 = 13, delta_0 = 8, num_steps = 200
  )

  set.seed(123)
  res <- do.call(gaussian_AD, default_parameters)
  set.seed(123)
  res2 <- do.call(gaussian_AD_2, default_parameters)

  for (i in names(res$d_beta)) {
    v1 <- colMeans_tail(res$d_beta[[i]])
    v2 <- colMeans_tail(res2$d_beta[[i]])
    expect_lt(sum(abs(v1 - v2)), 1e-10)
    print(rbind(implementation_1 = v1, implementation_2 = v2))
  }

  for (i in names(res$d_sigma2)) {
    v1 <- colMeans_tail(res$d_sigma2[[i]])
    v2 <- colMeans_tail(res2$d_sigma2[[i]])
    expect_lt(sum(abs(v1 - v2)), 1e-10)
    print(rbind(implementation_1 = v1, implementation_2 = v2))
  }
})
