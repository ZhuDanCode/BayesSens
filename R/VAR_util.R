# Create training data given a time series
train_data <- function(Y, lag) {
  nr <- nrow(Y)
  nc <- ncol(Y)
  X <- matrix(0, nr - lag, lag * nc)
  for (i in 1:lag) {
    X[, (1 + (i-1) * nc):(i * nc)] <- Y[(1 + lag - i):(nr - i), ]
  }
  list(X = cbind(1, X), Y = Y[(1 + lag):nr, ])
}
