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


# Collect objects by attributes
collect <- function(l0) {
  list_names <- names(l0[[1]])
  res <- vector("list", length(list_names))
  for (i in seq_along(list_names)) {
    res[[i]] <- purrr::map(l0, ~.x[[list_names[i]]]) %>% do.call(rbind, .)
  }
  purrr::set_names(res, list_names)
}
