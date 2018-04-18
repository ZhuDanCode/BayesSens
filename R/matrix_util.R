# Multiply a diagonal matrix by a matrix from the left without forming the
# diagonal matrix, M %*% diag(vec0).
left_multiply_D <- function(M, vec0) {
  t(t(M) * as.numeric(vec0))
}


neg_tx_otimes_x <- function(A) {
  -t(A) %x% A
}


kronecker_sp_1 <- function(A, B) {
  # kronecker_special_form_1: A %*% (I_n %x% B)
  n <- ncol(A) / nrow(B)
  C <- matrix(0, nrow(A), n * ncol(B))
  nr <- nrow(B)
  nc <- ncol(B)
  for (i in seq(n)) {
    C[, ((i-1)*nc + 1):(i*nc)] <- A[, ((i-1)*nr + 1):(i*nr), drop = F] %*% B
  }
  C
}


kronecker_sp_2 <- function(B, A) {
  # kronecker_special_form_1: (I_n %x% B) %*% A
  nc <- ncol(B)
  nr <- nrow(B)
  n <- nrow(A) / ncol(B)
  C <- matrix(0, nrow(B) * n, ncol(A))
  for (i in 1:n) {
    C[(1+(i-1)*nr):(i*nr), ] <- B %*% A[(1+(i-1)*nc):(i*nc), , drop = F]
  }
  C
}


# kronecker_sp_3b <- function(b, A) {
#   # kronecker_special_form_1: (b %x% I_n) %*% A
#   step <- nrow(A) / ncol(b)
#   assertthat::assert_that(step %% 1 == 0)
#   C <- matrix(0, nrow(b) * step, ncol(A))
#   for (i in seq(nrow(b))) {
#     for (j in seq(ncol(A))) {
#       for (k in seq(ncol(b))) {
#         C[(1 + step * (i-1)):(step * i), j] <-
#           C[(1 + step * (i-1)):(step * i), j] + b[i, k] * A[seq(step) + (k-1) * step, j]
#       }
#     }
#   }
#   C
# }


commutation_matrix <- function(r, c) {
  entries <- expand.grid(1:r, 1:c)
  src <- r * (entries[,2] - 1) + entries[,1]
  tgt <- c * (entries[,1] - 1) + entries[,2]
  Matrix::sparseMatrix(tgt, src, x = 1)
}


elimination_matrix <- function(n) {
  entries <- expand.grid(1:n, 1:n) %>%
    cbind(src = 1:(n^2)) %>%
    dplyr::filter(Var1 >= Var2) %>%
    cbind(tgt = 1:(0.5*n*(n+1)))
  Matrix::sparseMatrix(entries[,"tgt"], entries[,"src"], x = 1)
}
