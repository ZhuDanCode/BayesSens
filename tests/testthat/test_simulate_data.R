context("Testing simulation functions work for Gaussian, Student-t and SUR2.")

test_that("gaussian_data", {
  data0 <- gaussian_data(100, 10)
  expect_equal(length(data0), 4)
  expect_equal(names(data0), c("X", "y", "beta", "sigma"))
  expect_equal(any(is.na(data0$X)), FALSE)
  expect_equal(any(is.na(data0$y)), FALSE)
  expect_equal(any(is.na(data0$beta)), FALSE)
  expect_gt(data0$sigma, 0)
})

test_that("student_t_data", {
  data0 <- student_t_data(100, 10)
  expect_equal(length(data0), 4)
  expect_equal(names(data0), c("X", "y", "beta", "df"))
  expect_equal(any(is.na(data0$X)), FALSE)
  expect_equal(any(is.na(data0$y)), FALSE)
  expect_equal(any(is.na(data0$beta)), FALSE)
  expect_gt(data0$df, 0)
})

test_that("SUR2_data", {
  data0 <- SUR2_data(100, 10, 3)
  expect_equal(length(data0), 7)
  expect_equal(names(data0), c("Xy", "y", "beta", "Xs", "s", "gamma", "Sigma"))
  expect_equal(any(is.na(data0$Xy)), FALSE)
  expect_equal(any(is.na(data0$y)), FALSE)
  expect_equal(any(is.na(data0$beta)), FALSE)
  expect_equal(any(is.na(data0$Xs)), FALSE)
  expect_equal(any(is.na(data0$s)), FALSE)
  expect_equal(any(is.na(data0$gamma)), FALSE)
  expect_equal(any(is.na(data0$Sigma)), FALSE)
})
