# Optimised matrix multiplications
C_times_I_x_B <- function(C, B) {
  # kronecker_special_form_1: C %*% (I_n %x% B)
  n <- ncol(C) / nrow(B)
  Z <- matrix(0, nrow(C), n * ncol(B))
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
  Z <- matrix(0, n * nrow(B), ncol(C))
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
  Z <- matrix(0, nrow(A) * n, ncol(C))
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
  Z <- matrix(0, nrow(C), ncol(A) * n)
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
  s <- 1 + seq(0, by = n, length.out = n*q) %% (n*q - 1)
  s[n*q] <- n * q
  A[s,]
}
