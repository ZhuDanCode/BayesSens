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


# Multiply a diagonal matrix by a matrix from the left without forming the
# diagonal matrix, M %*% diag(vec0).
left_multiply_D <- function(M, vec0) {
  t(t(M) * vec0)
}
