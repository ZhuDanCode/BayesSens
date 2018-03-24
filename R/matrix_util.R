neg_tx_otimes_x <- function(A) {
  -t(A) %x% A
}


kronecker_sp_1 <- function(A, B) {
  # kronecker_special_form_1: A %*% (I_n %x% B)
  step <- nrow(B)
  do.call(cbind, purrr::map(
    seq(1, ncol(A), step),
    ~A[, .x:(.x + step - 1), drop = F] %*% B
  ))
}


kronecker_sp_2 <- function(B, A) {
  # kronecker_special_form_1: (I_n %x% B) %*% A
  step <- ncol(B)
  do.call(rbind, purrr::map(
    seq(1, nrow(A), step),
    ~B %*% A[.x:(.x + step - 1), , drop = F]
  ))
}


kronecker_sp_3 <- function(B, A) {
  # kronecker_special_form_1: (B %x% I_n) %*% A
  step <- nrow(A) / ncol(B)
  assertthat::assert_that(step %% 1 == 0)
  vec_by_rowsblock <- function(i) {
    s <- seq(1, nrow(A), step)
    purrr::map(
      1:ncol(B),
      ~B[i, .x] * A[s[.x]:(s[.x] + step - 1), , drop = F]
    ) %>%
      purrr::reduce(`+`)
  }
  1:nrow(B) %>%
    purrr::map(vec_by_rowsblock) %>%
    do.call(rbind, .)
}
