neg_tx_otimes_x <- function(A) {
  -t(A) %x% A
}


kronecker_sp_1 <- function(A, B) {
  # kronecker_special_form_1: A %*% (I_n %x% B)
  step <- ncol(B)
  do.call(cbind, purrr::map(
    seq(1, ncol(A), step),
    ~A[, .x:(.x + step - 1)] %*% B
  ))
}
