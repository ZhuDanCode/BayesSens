#' #' Compute the posterior statistics of the parameters
#' #' @param model0 Output from 'VAR_Gibbs'.
#' #' @param fun function; statistical function, default set to be 'mean'.
#' #' @export
#' VAR_posterior_stat <- function(model0, fun = mean) {
#'   n <- sqrt(ncol(model0$Sigma))
#'   res <- VAR_vec_to_model(apply(model0$beta, 2 , fun), n)
#'   res$Sigma <- matrix(apply(model0$Sigma, 2, fun), n, n)
#'   res
#' }
#'
#'
#' #' Compute the posterior forecast
#' #' @description This function takes each posterior sample and perform a n-step
#' #' forecast on the given time series.
#' #' @param model0 Output from 'VAR_Gibbs'.
#' #' @param data0 Matrix; each row is one observation, each column is one measurement / predictor.
#' #' @param n_step_ahead integer; number of future forecast to generate.
#' #' @return A array ARR where ARR[, , k] contains the forecast of the posterior sample at
#' #' iteration k. ARR[, , k] is of dimension D x n_step_ahead, where D is the dimension
#' #' of the multivariate time series.
#' #' @export
#' VAR_posterior_forecast <- function(model0, data0, n_step_ahead) {
#'   n <- ncol(data0)
#'   res <- 1:nrow(model0$beta) %>%
#'     purrr::map(~VAR_vec_to_model(model0$beta[.x, ], n)) %>%
#'     purrr::map(~tail(
#'       predict_VAR(n_step_ahead, data0, .x$b_0, .x$B),
#'       n_step_ahead
#'     ))
#'   array(unlist(res), dim = c(nrow(res[[1]]), ncol(res[[1]]), length(res)))
#' }
#'
#'
#' #' Summarise forecast posterior density
#' #' @param arr_forecast Output from VAR_posterior forecast
#' #' @param method function; summary statistical function, e.g. mean, quantile, etc.
#' #' @param ... Other parameters to be passed to the 'method'.
#' #' @return A list of summary statistics. The list should have length D, where D is
#' #' the dimension of the multivariate time series. Each component contains a matrix
#' #' with dimension K x n_step_ahead, where K is the dimension of the summary statistics
#' #' and n_step_ahead is the parameter set in "VAR_posterior_forecast".
#' #' @export
#' VAR_summarise_forecast <- function(arr_forecast, method, ...) {
#'   tmp_fun <- function(x) {
#'     apply(arr_forecast[x,,], 1, method, ...)
#'   }
#'   purrr::map(1:nrow(arr_forecast), tmp_fun)
#' }
