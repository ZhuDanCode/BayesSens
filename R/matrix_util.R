neg_tx_otimes_x <- function(A) {
  -t(A) %x% A
}

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

# commutation_matrix <- function(r, c) {
#   m0 <- matrix(0, r*c, r*c)
#   for (i in 1:r) {
#     for (j in 1:c) {
#       src <- r * (j - 1) + i
#       tgt <- c * (i - 1) + j
#       m0[tgt, src] <- 1
#     }
#   }
#   as(m0, "dgCMatrix")
# }

# elimination_matrix <- function(n) {
#   m0 <- matrix(0, 0.5 * n * (n + 1), n^2)
#   src <- 1
#   tgt <- 1
#   for (j in 1:n) {
#     for (i in 1:n) {
#       if (i >= j) {
#         m0[tgt, src] <- 1
#         tgt <- tgt + 1
#       }
#       src <- src + 1
#     }
#   }
#   as(m0, "dgCMatrix")
# }
