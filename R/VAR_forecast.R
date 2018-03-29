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
VAR_vec_to_model_coeff <- function(vec0, n) {
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


#' Compute the posterior statistics of the parameters
#' @param model0 Output from 'VAR_Gibbs'.
#' @param fun function; statistical function, default set to be 'mean'.
#' @export
VAR_posterior_stat <- function(model0, fun = mean) {
  n <- sqrt(ncol(model0$Sigma))
  res <- VAR_vec_to_model_coeff(apply(model0$beta, 2 , fun), n)
  res$Sigma <- matrix(apply(model0$Sigma, 2, fun), n, n)
  res
}


#' Compute the posterior statistics of the parameters
#' @param n_step_ahead integer; number of future forecast to generate.
#' @param data0 numeric matrix; the multivariate time series data.
#' The matrix must have at least as many columns as the lag, and the first
#' column should be the state at the first time point.
#' @param model0 Output from 'VAR_Gibbs'.
#' @export
VAR_posterior_forecast <- function(n_step_ahead, data0, model0) {
  n <- nrow(data0)
  last_n_cols <- function(m0, n) { m0 %>% t() %>% tail(n) %>% t() }
  res <- 1:nrow(model0$beta) %>%
    purrr::map(~VAR_vec_to_model_coeff(model0$beta[.x, ], n)) %>%
    purrr::map(~last_n_cols(
      predict_VAR(n_step_ahead, data0, .x$b_0, .x$B),
      n_step_ahead
    ))
  array(unlist(res), dim = c(nrow(res[[1]]), ncol(res[[1]]), length(res)))
}


#' Predict future values given a VAR model
#' @param n_step_ahead integer; number of future forecast to generate.
#' @param data0 numeric matrix; the multivariate time series data.
#' The matrix must have at least as many columns as the lag, and the first
#' column should be the state at the first time point.
#' @param b_0 numeric vector; the intercept.
#' @param B list of numeric matrices; the regression coefficients.
#' @export
predict_VAR <- function(n_step_ahead, data0, b_0, B) {
  nr <- nrow(data0)
  nc <- ncol(data0)
  data1 <- cbind(data0, matrix(0, nr, n_step_ahead))
  lag <- length(B)
  for (i in nc + seq(n_step_ahead)) {
    data1[, i] <- fit(b_0, B, data1[, (i-lag):(i-1)])
  }
  data1
}


#' Generate fitted values given a VAR model
#' @param data0 numeric matrix; the multivariate time series data.
#' The matrix must have at least as many columns as the lag, and the first
#' column should be the state at the first time point.
#' @param b_0 numeric vector; the intercept.
#' @param B list of numeric matrices; the regression coefficients.
#' @export
fitted_VAR <- function(data0, b_0, B) {
  data1 <- data0
  lag <- length(B)
  for (i in (lag + 1):ncol(data0)) {
    data1[,i] <- fit(b_0, B, data0[,(i-lag):(i-1)])
  }
  data1
}


fit <- function(b_0, B, y_t) {
  z <- b_0
  for (i in seq_along(B)) {
    z <- z + B[[i]] %*% y_t[,i]
  }
  z
}


#' Summarise forecast posterior density
#' @param arr_forecast Output from VAR_posterior forecast
#' @param method function; summary statistical function, e.g. mean, quantile, etc.
#' @param ... Other parameters to be passed to the 'method'.
#' @export
summarise_forecast <- function(arr_forecast, method, ...) {
  tmp_fun <- function(x) {
    apply(arr_forecast[x,,], 1, method, ...)
  }
  purrr::map(1:nrow(arr_forecast), tmp_fun)
}
