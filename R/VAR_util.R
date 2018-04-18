# VAR util functions
#' Vectorise all the coefficients appeared in the model
#' @param model0 A VAR model; the second component of the output from 'simulate_VAR'.
#' @export
VAR_model_to_vec <- function(model0) {
  matrixcalc::vec(t(cbind(model0$b_0, do.call(cbind, model0$B))))
}


#' Convert vector back to model coefficients
#' @param vec0 Numeric vector; vectorised model coefficients.
#' @param n integer; dimension of the data.
#' @export
VAR_vec_to_model <- function(vec0, n) {
  m0 <- matrix(vec0, ncol = n, byrow = F)
  s <- seq(2, nrow(m0), n)
  list(
    b_0 = m0[1,],
    B = purrr::map2(s, s+n-1, ~t(m0[.x:.y, , drop = F]))
  )
}


#' Convert the data into a format ready for modelling
#' @param data0 Matrix; the number of columns is equal to the number of time points,
#' the number of rows is equal to the dimension of the observation. One can get
#' such data from the first component of the output from 'simulate_VAR'.
#' @param p integer; time series lag.
#' @export
VAR_model_matrix <- function(data0, p) {
  eff_T <- ncol(data0) - p
  n <- nrow(data0)
  I_n <- diag(n)
  purrr::map2(
    .x = seq(eff_T), .y = p + seq(eff_T) - 1,
    .f = ~ I_n %x% t( c(1, matrixcalc::vec(data0[, .y:.x, drop = F])) )
  ) %>%
    do.call(rbind, .)
}


# Create training data given a time series
#' @keywords internal
train_data <- function(Y, lag) {
  nr <- nrow(Y)
  nc <- ncol(Y)
  X <- matrix(0, nr - lag, lag * nc)
  for (i in 1:lag) {
    X[, (1 + (i-1) * nc):(i * nc)] <- Y[(1 + lag - i):(nr - i), ]
  }
  list(X = cbind(1, X), Y = Y[(1 + lag):nr, ])
}
