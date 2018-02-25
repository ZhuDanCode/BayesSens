#' Infinitesimal Perturbation Analysis
#' @param model0 Function; any one of *-Gibbs functions.
#' @param X Numeric matrix; the covariates.
#' @param y Numeric vector; the response variable.
#' @param param List; a list of named parameters.
#' @param change String; name of the parameter to change.
#' @param pos Integer; if the parameter is a vector, then provide the position
#' of the specific parameter.
#' @param h Number; the step size for the sensitivitiy analysis.
#' @param stat_fun Function; the summary function of the posterior distribution.
#' @examples
#' \dontrun{
#' b_0 <- rnorm(p)
#' B_0 <- pdmatrix(p)$Sigma
#' alpha_0 <- 13
#' delta_0 <- 8
#' # Analyse the sensitivity of posterior mean w.r.t. parameter 'alpha_0'.
#' h <- 0.01
#' stat_fun <- mean
#' IPA(
#'   gaussian_Gibbs, data0$X, data0$y,
#'   param = list(b_0 = b_0, B_0 = B_0, alpha_0 = alpha_0, delta_0 = delta_0),
#'   change = "alpha_0", h = h, stat_fun = stat_fun
#' )
#' }
#' @export
#'
# IPA <- function(model0, X, y, param, change, pos = 1, h, stat_fun = mean) {
#   new_param <- param
#   new_param[[change]][pos] <- new_param[[change]][pos] + h
#   data0 <- list(X = X, y = y)
#
#   base_res <- do.call(model0, append(data0, param))
#   compare_res <- do.call(model0, append(data0, new_param))
#
#   compute_stat <- function(x) {
#     if (is.null(nrow(x))) return(stat_fun(x))
#     apply(x, 2, stat_fun)
#   }
#   compute_sense <- function(x, y) {
#     (compute_stat(x) - compute_stat(y)) / h
#   }
#
#   purrr::map2(base_res, compare_res, ~compute_sense(.x, .y)) %>%
#     setNames(paste0(
#       names(base_res), "_sensitivity_to_", change, ifelse(pos == 1, "", pos)
#     ))
# }
