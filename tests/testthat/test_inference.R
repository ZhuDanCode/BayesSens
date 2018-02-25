library(autoBayes)
test_that("gaussian_Gibbs", {
  skip_on_cran()
  n <- 1000
  p <- 5
  set.seed(123)
  data0 <- gaussian_data(n, p, intercept = TRUE)
  res <- gaussian_Gibbs(data0$X, data0$y,
    b_0 = rnorm(p+1), B_0 = pdmatrix(p+1)$Sigma,
    alpha_0 = 13, delta_0 = 8
  )
  # check with frequentist approach
  glm_model0 <- glm(data0$y ~ data0$X - 1, family = gaussian())
  glm_beta <- coefficients(glm_model0)         # beta
  glm_sigma <- sqrt(deviance(glm_model0) / n)  # sigma
  # test
  expect_lt(sum(abs(apply(res$beta, 2, mean) - glm_beta)), 1 / sqrt(n))
  expect_lt(abs(mean(res$sigma) - glm_sigma), 1 / sqrt(n))
})
