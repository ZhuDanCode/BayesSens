test_that("gaussian_data", {
  data0 <- gaussian_data(100, 10)
  expect_equal(length(data0), 4)
  expect_equal(names(data0), c("X", "y", "beta", "sigma"))
  expect_equal(any(is.na(data0$X)), FALSE)
  expect_equal(any(is.na(data0$y)), FALSE)
  expect_equal(any(is.na(data0$beta)), FALSE)
  expect_gt(data0$sigma, 0)
})
