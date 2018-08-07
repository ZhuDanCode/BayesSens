apply_chain <- function(expr_fun, hyperparameter) {
  if (!is.function(expr_fun)) return(0)
  res <- purrr::map(hyperparameter, expr_fun)
  names(res) <- hyperparameter
  res
}


d_chol <- function(L, dA, I_nn, K_nn, I_n, E_n, fac_1) {
  # LL^T = A
  n <- nrow(L)
  if (missing(I_n)) I_n <- Matrix::Diagonal(n)
  if (missing(I_nn)) I_nn <- Matrix::Diagonal(n^2)
  if (missing(K_nn)) K_nn <- commutation_matrix(n, n)
  if (missing(E_n)) E_n <- elimination_matrix(n)
  D_n <- Matrix::t(E_n)
  if (missing(fac_1)) fac_1 <- E_n + A_times_K_nq(E_n, n, n)  # E_n %*% (I_nn + K_nn)
  fac_2 <- as.matrix(
    D_n %*% solve(fac_1 %*% kronecker_sp_3_cpp(L, as.matrix(D_n))) %*% E_n)
  # apply_chain(. %>% {fac_2 %*% dA[[.]]})
  hyperparam <- names(dA)
  apply_chain(. %>% {eigenMapMatMult(fac_2, as.matrix(dA[[.]]))}, hyperparam)
}


d_transpose <- function(X, dX, K_nq) {
  apply_chain(. %>% {K_nq_times_A(dX[[.]], nrow(X), ncol(X))}, names(dX))
}


d_XXT <- function(X, dX, I_nn, K_nn, I_n) {
  n <- nrow(X)
  if (missing(I_n)) I_n <- Matrix::Diagonal(n)
  # if (missing(I_nn)) I_nn <- Matrix::Diagonal(n^2)
  # if (missing(K_nn)) K_nn <- commutation_matrix(n, n)
  hyperparam <- names(dX)
  if (n < 50) {
    m0 <- X %x% I_n
    f1 <- m0 + K_nq_times_A(m0, n, n)
    # f1 <- (I_nn + K_nn) %*% (X %x% I_n)
  } else {
    f1 <- C_times_A_x_I(I_nn + K_nn, X)
  }
  apply_chain(. %>% {f1 %*% dX[[.]]}, hyperparam)
}


d_kronecker <- function(A, dA, B, dB, fac_1, I_mn, I_pq, I_n, K_qm, I_p) {
  # apply_chain(. %>% {A %x% dB[[.]] + dA[[.]] %x% B})
  hyperparam <- names(dA)
  m <- nrow(A)
  n <- ncol(A)
  p <- nrow(B)
  q <- ncol(B)
  if (missing(fac_1)) {
    if (missing(I_n)) I_n <- Matrix::Diagonal(n)
    if (missing(K_qm)) K_qm <- commutation_matrix(q, m)
    if (missing(I_p)) I_p <- Matrix::Diagonal(p)
    if (missing(I_mn)) I_mn <- Matrix::Diagonal(m*n)
    if (missing(I_pq)) I_pq <- Matrix::Diagonal(p*q)
    fac_1 <- (I_n %x% K_qm %x% I_p)
  }
  ind <- case(dA, dB)
  expr_fun <- switch(ind,
     '0' = . %>% {fac_1 %*% (as.numeric(A) %x% dB[[.]] + dA[[.]] %x% as.numeric(B))},
     '1' = . %>% { fac_1 %*% (as.numeric(A) %x% dB[[.]]) },
     '2' = . %>% { fac_1 %*% (dA[[.]] %x% as.numeric(B)) },
     '3' = 0)
  hyperparam <- switch(ind, '2' = names(dA), '3' = 0, names(dB))
  apply_chain(expr_fun, hyperparam)
  # return(apply_chain(. %>% { fac_1 %*% ((I_mn %x% as.numeric(B)) %*% dA[[.]]) }))
  # return(apply_chain(. %>% { fac_1 %*% (I_x_B_times_A(as.numeric(B), dA[[.]])) }))
  # return(apply_chain(. %>% { fac_1 %*% ((as.numeric(A) %x% I_pq) %*% dB[[.]]) }))
  # return(apply_chain(. %>% { fac_1 %*% (I_x_B_times_A(as.numeric(A), dB[[.]])) }))
  # apply_chain(. %>% {fac_1 %*%
  # ( (as.numeric(A) %x% I_pq) %*% dB[[.]] + (I_mn %x% as.numeric(B)) %*% dA[[.]] )})
  # apply_chain(. %>% {fac_1 %*%
  # ( kronecker_sp_3(as.numeric(A), dB[[.]]) + I_x_B_times_A(as.numeric(B), dA[[.]]) )})
}


d_sum <- function(dA, dB) {
  if (is_zero(dA) && is_zero(dB)) {return(0)}
  if (is_zero(dB)) { return(dA) }
  if (is_zero(dA)) { return(dB) }
  apply_chain(. %>% {dA[[.]] + dB[[.]]}, names(dA))
}


d_product <- function(A, dA, B, dB, I_a, I_b) {
  if (is.vector(A)) A <- as.matrix(A)
  if (is.vector(B)) B <- as.matrix(B)
  if (missing(I_a)) I_a <- Matrix::Diagonal(nrow(A))
  if (missing(I_b)) I_b <- Matrix::Diagonal(ncol(B))
  ind <- case(dA, dB)
  # expr_fun <- switch(ind,
  #    '0' = . %>% {(I_b %x% A) %*% dB[[.]] + (t(B) %x% I_a) %*% dA[[.]]},
  #    '1' = . %>% {(I_b %x% A) %*% dB[[.]]},
  #    '2' = . %>% {(t(B) %x% I_a) %*% dA[[.]]},
  #    '3' = 0)
  expr_fun <- switch(ind,
     '0' = . %>% {I_x_B_times_C(A, dB[[.]]) + A_x_I_times_C(t(B), dA[[.]])},
     '1' = . %>% {I_x_B_times_C(A, dB[[.]])},
     '2' = . %>% {A_x_I_times_C(t(B), dA[[.]])},
     '3' = 0)
  hyperparam <- switch(ind, '2' = names(dA), '3' = 0, names(dB))
  apply_chain(expr_fun, hyperparam)
}


d_minus <- function(dA, dB) {
  ind <- case(dA, dB)
  expr_fun <- switch(ind,
     '0' = . %>% {dA[[.]] - dB[[.]]},
     '1' = . %>% {- dB[[.]]},
     '2' = . %>% {dA[[.]]},
     '3' = 0)
  hyperparam <- switch(ind, '2' = names(dA), '3' = 0, names(dB))
  apply_chain(expr_fun, hyperparam)
}


d_inv <- function(inv_A, dA, fac_1) {
  if (missing(fac_1)) fac_1 <- - t(inv_A) %x% inv_A
  # expr <- . %>% {fac_1 %*% dA[[.]]}
  expr <- . %>% {eigenMapMatMult(fac_1, as.matrix(dA[[.]]))}
  apply_chain(expr, names(dA))
}


d_Gamma <- function(g, alpha) {
  f <- function(t) { log(t) * dgamma(t, alpha, 1) }
  num_1 <- integrate(f, 0, g)$value
  num_2 <- digamma(alpha) * pgamma(g, alpha, 1)
  - (num_1 - num_2) / dgamma(g, alpha, 1)
}


d_constant_inv <- function(x, dx) {
  x <- 1 / x
  apply_chain(. %>% {-x^(-2) * dx[[.]]}, names(dx))
}


d_constant_multiply_matrix <- function(c, dc, X, dX) {
  ind <- case(dc, dX)
  expr_fun <- switch(ind,
    '0' = . %>% {as.numeric(X) %*% dc[[.]] + c * dX[[.]]},
    '1' = . %>% {c * dX[[.]]},
    '2' = . %>% {as.numeric(X) %*% dc[[.]]},
    '3' = 0)
  hyperparam <- switch(ind, '2' = names(dc), '3' = 0, names(dX))
  apply_chain(expr_fun, hyperparam)
}


d_constant_multiply_constant <- function(a, da, b, db) {
  ind <- case(da, db)
  expr_fun <- switch(ind,
    '0' = . %>% {a * db[[.]] + b * da[[.]]},
    '1' = . %>% {a * db[[.]]},
    '2' = . %>% {b * da[[.]]},
    '3' = 0)
  hyperparam <- switch(ind, '2' = names(da), '3' = 0, names(db))
  apply_chain(expr_fun, hyperparam)
}


d_constant_divide_constant <- function(a, da, b, db) {
  d_constant_multiply_constant(
    a, da, 1 / b, d_constant_inv(1 / b, db)
  )
}


is_zero <- function(x) {
  !is.list(x) && (x == 0)
}


case <- function(x, y) {
  # 0: list & list, 1: 0 & list, 2: list & 0, 3: 0 & 0
  as.character(ifelse(is_zero(x), 1, 0) + ifelse(is_zero(y), 2, 0))
}
