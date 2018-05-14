neg_tx_otimes_x <- function(A) {
  -t(A) %x% A
}


commutation_matrix <- function(r, c) {
  m0 <- matrix(0, r*c, r*c)
  for (i in 1:r) {
    for (j in 1:c) {
      src <- r * (j - 1) + i
      tgt <- c * (i - 1) + j
      m0[tgt, src] <- 1
    }
  }
  as(m0, "dgCMatrix")
}


elimination_matrix <- function(n) {
  m0 <- matrix(0, 0.5 * n * (n + 1), n^2)
  src <- 1
  tgt <- 1
  for (j in 1:n) {
    for (i in 1:n) {
      if (i >= j) {
        m0[tgt, src] <- 1
        tgt <- tgt + 1
      }
      src <- src + 1
    }
  }
  as(m0, "dgCMatrix")
}
