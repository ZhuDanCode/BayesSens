push <- function(l0, item) {
  l0[[length(l0) + 1]] <- item
  l0
}


#' This is a function to generate covariance matrix.
#' @details Alias for 'clusterGeneration::genPositiveDefMat'
#' @param ... Parameters to be passed to 'genPositiveDefMat'.
#' See '?clusterGeneration::genPositiveDefMat' for documentation.
#' @examples
#' # Create a random 4 x 4 covariance matrix
#' cov_mat <- pdmatrix(dim = 4)$Sigma
#' print(cov_mat)
#' @export
pdmatrix <- function(...) {
  clusterGeneration::genPositiveDefMat(...)
}


zeros <- function(nr = 1, nc = 1) {
  matrix(numeric(nr * nc), nr, nc)
}


init_differential <- function(len0, l_b0) {
  list(
    d_b0 = zeros(len0, l_b0),
    d_B0 = zeros(len0, l_b0^2),
    d_alpha0 = zeros(len0, 1),
    d_delta0 = zeros(len0, 1),
    d_sigma2_0 = zeros(len0, 1)
  )
}
