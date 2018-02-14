#' Plot the posterior distribution of the parameters
#' @param res Output from any *-Gibbs functions.
#' @export
plot_posterior <- function(res) {
  len <- num_param(res)
  nr <- floor(sqrt(len))
  nc <- ceiling(len / nr)
  par(mfrow = c(nr, nc))

  plot_hist <- function(x, lab) {
    hist(x, 20, probability = T, xlab = lab, main = 'Posterior distribution')
  }
  is_single_col <- . %>% {is.null(ncol(.))}

  for (i in seq_along(res)) {
    param <- res[[i]]
    if (is_single_col(param)) {
      plot_hist(param, names(res)[i])
    } else {
      for (j in 1:ncol(param)) {
        plot_hist(param[,j], paste0(names(res)[i], '_', j))
      }
    }
  }
  par(mfrow = c(1,1))
}


# Number of parameters in the output from *-Gibbs functions
num_param <- function(res) {
  res %>% purrr::map_dbl(~max(ncol(.x), 1)) %>% sum()
}
