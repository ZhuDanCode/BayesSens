#' Predict future values given a VAR model
#' @param n_step_ahead integer; number of future forecast to generate.
#' @param data0 Matrix; each row is one observation, each column is one measurement / predictor.
#' @param b_0 numeric vector; the intercept.
#' @param B list of numeric matrices; the regression coefficients.
#' @export
predict_VAR <- function(n_step_ahead, data0, b_0, B) {
  data0 <- t(data0)
  data1 <- cbind(data0, matrix(0, nrow(data0), n_step_ahead))
  lag <- length(B)
  for (i in ncol(data0) + seq(n_step_ahead)) {
    data1[, i] <- fit(b_0, B, data1[, (i-1):(i-lag)])
  }
  t(data1)
}


#' Generate fitted values given a VAR model
#' @param data0 Matrix; each row is one observation, each column is one measurement / predictor.
#' @param b_0 numeric vector; the intercept.
#' @param B list of numeric matrices; the regression coefficients.
#' @export
fitted_VAR <- function(data0, b_0, B) {
  data0 <- t(data0)
  data1 <- data0
  lag <- length(B)
  for (i in (lag + 1):ncol(data0)) {
    data1[,i] <- fit(b_0, B, data0[, (i-1):(i-lag)])
  }
  t(data1)
}


fit <- function(b_0, B, y_t) {
  z <- b_0
  for (i in seq_along(B)) {
    z <- z + B[[i]] %*% y_t[,i]
  }
  z
}
