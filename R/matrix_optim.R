neg_tx_otimes_x <- function(A) {
  -t(A) %x% A
}

commutation_matrix <- function(r, c) {
  if (missing(c)) c <- r
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

# Optimised matrix multiplications
C_times_I_x_B <- function(C, B) {
  # kronecker_special_form_1: C %*% (I_n %x% B)
  n <- ncol(C) / nrow(B)
  Z <- memo_zero_matrix(nrow(C), n * ncol(B))
  nr <- nrow(B)
  nc <- ncol(B)
  for (i in seq(n)) {
    Z[, ((i-1)*nc + 1):(i*nc)] <- C[, (1 + (i-1)*nr):(i*nr), drop = F] %*% B
  }
  Z
}

I_x_B_times_C <- function(B, C) {
  # kronecker_special_form_2: (I_n %x% B) %*% C
  n <- nrow(C) / ncol(B)
  Z <- memo_zero_matrix(n * nrow(B), ncol(C))
  nc <- ncol(B)
  nr <- nrow(B)
  for (i in 1:n) {
    Z[(1+(i-1)*nr):(i*nr), ] <- as.matrix(B %*% C[(1+(i-1)*nc):(i*nc), , drop = F])
  }
  Z
}

A_x_I_times_C <- function(A, C) {
  # kronecker_special_form_3: (A %x% I_n) %*% C
  n <- nrow(C) / ncol(A)
  Z <- memo_zero_matrix(nrow(A) * n, ncol(C))
  nc <- ncol(A)
  nr <- nrow(A)
  for (i in 1:nrow(A)) {
    for (j in 1:ncol(A)) {
      Z[(1+(i-1)*n):(i*n), ] <- Z[(1+(i-1)*n):(i*n), ] +
        diag_v0_times_A(rep(A[i,j], n), C[(1+(j-1)*n):(j*n), , drop = F])
    }
  }
  Z
}

C_times_A_x_I <- function(C, A) {
  # kronecker_special_form_4: C %*% (A %x% I_n)
  n <- ncol(C) / nrow(A)
  Z <- memo_zero_matrix(nrow(C), ncol(A) * n)
  nc <- ncol(A)
  nr <- nrow(A)
  for (i in 1:nrow(A)) {
    for (j in 1:ncol(A)) {
      Z[, (1+(j-1)*n):(j*n)] <- Z[, (1+(j-1)*n):(j*n)] +
        A_times_diag_v0(as.matrix(C[, (1+(i-1)*n):(i*n), drop = F]), rep(A[i,j], n))
    }
  }
  Z
}

A_times_diag_v0 <- function(A, v0) {
  t(t(A) * v0)
}

diag_v0_times_A <- function(v0, A) {
  v0 * A
}

A_times_K_nq <- function(A, n, q) {
  s <- 1 + seq(0, by = q, length.out = n*q) %% (n*q - 1)
  s[n*q] <- n * q
  A[,s]
}

K_nq_times_A <- function(A, n, q) {
  n <- as.numeric(n)
  q <- as.numeric(q)
  s <- 1 + seq(0, by = n, length.out = n*q) %% (n*q - 1)
  s[n*q] <- n * q
  A[s,]
}


memoize <- function(f) {
  table0 <- list()
  char_hash <- function(...) {
    paste0(as.character(list(...)), collapse = ",")
  }
  function(...) {
    char_x <- char_hash(...)
    res <- table0[[char_x]]
    if (is.null(res)) {
      res <- f(...)
      table0[[char_x]] <<- res
    }
    return(res)
  }
}


zero_matrix <- function(nr, nc) { matrix(0, nr, nc) }
memo_zero_matrix <- memoize(zero_matrix)
memo_Diagonal <- memoize(Matrix::Diagonal)
