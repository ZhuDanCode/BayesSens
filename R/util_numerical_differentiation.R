#' Numerical differentiation
#' @param f a function to be differentiated numerically.
#' @param fixed_param a list of fixed parameters to be called upon f.
#' @param vary_param a list of parameters of which the deriatives are seeked.
#' @param sum_stat a function to compute summary statistics.
#' @param delta a positive number, the amount of perturbation.
num_diff <- function(f, fixed_param, vary_param, sum_stat, delta, seed) {
  if (!missing(seed)) set.seed(seed)
  default_param <- append(fixed_param, vary_param)
  base_case <- do.call(f, default_param)

  vec_param <- slice_list(vary_param)
  res <- vector("list", length(vec_param))
  for (i in seq_along(vec_param)) {
    new_param <- vec_param
    new_param[i] <- new_param[i] + delta
    new_param <- stack_list(new_param)
    new_param <- append(fixed_param, new_param)

    if (!missing(seed)) set.seed(seed)
    new_case <- do.call(f, new_param)
    res[[i]] <- (sum_stat(new_case) - sum_stat(base_case)) / delta
    print(sprintf("Complete %d / %d", i, length(vec_param)))
  }
  res
}

