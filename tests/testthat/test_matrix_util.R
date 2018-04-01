testthat::test_that("test_matrix_special_form_1",{
  for (i in 1:30) {
    m <- sample(3:10, 1)
    d <- sample(3:10, 1)
    p <- sample(3:10, 1)
    q <- sample(1:3, 1)
    A <- matrix(rnorm(m * d * p), m, d * p)
    B <- matrix(rnorm(p * q), p, q)
    I_d <- diag(d)
    testthat::expect_equal(
      A %*% (I_d %x% B),
      kronecker_sp_1(A, B)
    )
  }
})


testthat::test_that("test_matrix_special_form_2",{
  for (i in 1:30) {
    m <- sample(3:10, 1)
    d <- sample(3:10, 1)
    p <- sample(3:10, 1)
    q <- sample(1:3, 1)
    A <- matrix(rnorm(m * d * q), d * q, m)
    B <- matrix(rnorm(p * q), p, q)
    I_d <- diag(d)
    testthat::expect_equal(
      (I_d %x% B) %*% A,
      kronecker_sp_2(B, A)
    )
  }
})


testthat::test_that("test_matrix_special_form_3",{
  for (i in 1:30) {
    m <- sample(3:10, 1)
    d <- sample(3:10, 1)
    p <- sample(3:10, 1)
    q <- sample(1:3, 1)
    A <- matrix(rnorm(m * d * q), d * q, m)
    B <- matrix(rnorm(p * q), p, q)
    I_d <- diag(d)
    testthat::expect_equal(
      (B %x% I_d) %*% A,
      kronecker_sp_3_cpp(B, A)
    )
  }
})


testthat::test_that("Commutation Matrix", {
  for (i in 1:10) {
    m <- sample(2:10, 1)
    n <- sample(2:10, 1)
    testthat::expect_equal(
      commutation_matrix(m,n),
      as(matrixcalc::commutation.matrix(m, n), "dgCMatrix")
    )
  }
})

test_that("Elimination matrix", {
  for (i in 1:10) {
    m <- sample(2:15, 1)
    expect_equal(
      elimination_matrix(m),
      as(matrixcalc::elimination.matrix(m), "dgCMatrix")
    )
  }
})
