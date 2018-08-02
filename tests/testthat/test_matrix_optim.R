context("Testing optimised matrix operations produce correct results")

testthat::test_that("Left and Right Multiplication by a diagonal matrix", {
  n <- 50
  A <- matrix(rnorm(n^2), n, n)
  v0 <- rnorm(n)

  testthat::expect_equal(
    A_times_diag_v0(A, v0),
    A %*% diag(v0)
  )

  testthat::expect_equal(
    diag_v0_times_A(v0, A),
    diag(v0) %*% A
  )
})

testthat::test_that("Kronecker product with diagonal matrices", {
  # C %*% (I_n %x% B)
  n <- sample(1:50, 100, replace = TRUE)
  k <- sample(1:10, 100, replace = TRUE)
  p <- n * k
  for(i in 1:100) {
    C <- matrix(rnorm(p[i]^2), p[i], p[i])
    B <- matrix(rnorm(k[i]^2), k[i], k[i])
    I_n <- diag(n[i])
    testthat::expect_equal(
      C_times_I_x_B(C, B),
      C %*% (I_n %x% B)
    )
  }

  # (I_n %x% B) %*% C
  n <- sample(1:50, 100, replace = TRUE)
  k <- sample(1:10, 100, replace = TRUE)
  p <- n * k
  for(i in 1:100) {
    I_n <- diag(n[i])
    B <- matrix(rnorm(k[i]^2), k[i], k[i])
    C <- matrix(rnorm(p[i]^2), p[i], p[i])
    testthat::expect_equal(
      I_x_B_times_C(B, C),
      (I_n %x% B) %*% C
    )
  }

  # (A %x% I_n) %*% C
  n <- sample(1:50, 100, replace = TRUE)
  k <- sample(1:10, 100, replace = TRUE)
  p <- n * k
  for(i in 1:100) {
    A <- matrix(rnorm(k[i]^2), k[i], k[i])
    I_n <- diag(n[i])
    C <- matrix(rnorm(p[i]^2), p[i], p[i])
    testthat::expect_equal(
      A_x_I_times_C(A, C),
      (A %x% I_n) %*% C
    )
  }

  # C %*% (A %x% I_n)
  n <- sample(1:50, 100, replace = TRUE)
  k <- sample(1:10, 100, replace = TRUE)
  p <- n * k
  for(i in 1:100) {
    C <- matrix(rnorm(p[i]^2), p[i], p[i])
    A <- matrix(rnorm(k[i]^2), k[i], k[i])
    I_n <- diag(n[i])
    testthat::expect_equal(
      C_times_A_x_I(C, A),
      C %*% (A %x% I_n)
    )
  }
})
