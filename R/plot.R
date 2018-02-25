#' Plot the posterior distribution of the parameters
#' @param res Output from any *-Gibbs functions.
#' @export
plot_posterior <- function(res) {
  plot_hist <- function(x, lab) {
    hist(x, 20, probability = T, xlab = lab, main = 'Posterior distribution')
  }
  is_single_col <- . %>% {is.null(ncol(.))}

  reshape_subplots(num_param(res))
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
  reshape_subplots()
}


# Number of parameters in the output from *-Gibbs functions
num_param <- function(res) {
  res %>% purrr::map_dbl(~max(ncol(.x), 1)) %>% sum()
}


reshape_subplots <- function(len) {
  if (missing(len)) {
    par(mfrow = c(1,1))
  } else {
    nr <- floor(sqrt(len))
    nc <- ceiling(len / nr)
    par(mfrow = c(nr, nc))
  }
  invisible(NULL)
}


# plot_col <- function(df0) {
#   nc <- ncol(df0)
#   reshape_subplots(min(6, nc))
#   for (i in 1:nc) {
#     plot(df0[,i], type = 'l')
#     if (i %% 6 == 0) {
#       input0 <- readline("Continue (Y or N)? ")
#       if (input0 == "N") {
#         return(reshape_subplots())
#       }
#     }
#   }
#   reshape_subplots()
# }
