# context("Testing conversion between different data structures")
#
# testthat::test_that("vec <-> matrix conversion", {
#   x <- vec_to_mat(1:10, c(2, 5))
#   expect_equal(x, vec_to_mat(mat_to_vec(x)))
#   x <- vec_to_arr(1:18, c(2, 2, 3))
#   expect_equal(x, vec_to_arr(arr_to_vec(x)))
# })
#
# testthat::test_that("vec <-> list <vec> conversion", {
#   x <- list(1:3, 1:8)
#   expect_equal(x, vec_to_listv(listv_to_vec(x)))
# })
#
# testthat::test_that("vec <-> list <array> conversion", {
#   x <- list(1:3, matrix(1:4, 2, 2))
#   expect_equal(x, stack_list(slice_list(x)))
# })
