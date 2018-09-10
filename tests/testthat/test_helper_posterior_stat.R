testthat::test_that("Test Jacobian and Cumulative Stat", {
  x <- rnorm(100*10*3)
  dim(x) <- c(100, 10, 3)
  testthat::expect_equal(
    jacobian(x, cumulative_stat, fun = mean)[,1,1],
    cumulative_stat(x[,1,1], mean)
  )
})
