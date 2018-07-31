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


#' This is a function to generate matrix with eigenvalue less than 1 based on
#' some heuristic.
#' @param dim integer; dimension of the matrix.
#' @examples
#' # Create a random 4 x 4 matrix with eigenvalue less than 1.
#' mat <- s_eig_matrix(4)
#' print(mat)
#' eigen(mat)$value
#' @export
s_eig_matrix <- function(dim) {
  matrix(rnorm(dim^2), dim, dim) %*% diag(runif(dim, 0, 1/dim)) %*% matrix(rnorm(dim^2), dim, dim)
}


# Create zero matrix
zeros <- function(nr = 1, nc = 1) {
  matrix(numeric(nr * nc), nr, nc)
}


# Initialise internal input for Autodiff function
init_differential <- function(len0, vec0, names0) {
  # initialise list of zero matrices
  # len0 should be the dimension of the numerator
  # vec0 should be the vector of dimensions of the denominator (parameters)
  # names0 is the vector of names of the denominator (parameters)
  names0 <- paste0("d_", names0)
  vec0 %>%
    purrr::map(~zeros(len0, .x)) %>%
    set_names(names0)
}


# Tidy internal output for Autodiff function
# input: dlist0 contains N lists of matrices named by vec0.
# output: a named list of (N-row) matrices.
tidy_list <- function(dlist0) {
  extract_rbind <- function(attr0) {
    dlist0 %>%
      purrr::map(~t(as.numeric(.x[[attr0]]))) %>%
      do.call(rbind, .)
  }
  vec0 <- names(dlist0[[1]])
  vec0 %>%
    purrr::map(extract_rbind) %>%
    purrr::set_names(vec0)
}
