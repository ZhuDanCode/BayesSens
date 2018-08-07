# This file contains helper functions for testing

# Summary statistics - Posterior mean
colMeans_tail <- function(X, p = 0.9) {
  colMeans(tail(X, max(round(nrow(X) * p), 1)))
}

mean_tail <- function(v0, p = 0.9) {
  mean(tail(v0, max(round(length(v0) * p), 1)))
}


# Posterior mean calculation for numerical differentiation and auto-differentiation
num_diff <- function(x, y, h, ...) {
  # x, y: matrices with 'num_steps' rows and 'dim(parameters)' cols.
  tail_mean_fun <- if_else(
    test = (is.matrix(x) && is.matrix(y)),
    yes = colMeans_tail,  # matrix case
    no = mean_tail        # vector case
  )
  (tail_mean_fun(x, ...) - tail_mean_fun(y, ...)) / h
}

auto_diff <- function(res, num, dem) {
  apply(get_sensitivity(res, num, dem, T), 2:3, mean_tail)
}
