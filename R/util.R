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


# Create zero matrix
zeros <- function(nr = 1, nc = 1) {
  matrix(numeric(nr * nc), nr, nc)
}


#' Initialise internal input for Autodiff function
#' @description This function initialises list of zero matrices.
#' @param len0 the dimension of the numerator
#' @param vec0 the vector of dimensions of the denominator (parameters)
#' @param names0 the vector of names of the denominator (parameters)
init_differential <- function(len0, vec0, names0) {
  map_named(vec0, ~zeros(len0, .x), paste0("d_", names0))
}

differential_matrix <- function(n) {
  res <- commutation_matrix(n, n)
  for (i in seq(n^2)) {
    res[i,i] <- 1
  }
  res
}

#' Tidy internal output for Autodiff function
#' @param dlist0 expects (a list of) N lists of matrices named by vec0.
#' @return a named list of (N-row) matrices.
tidy_list <- function(dlist0) {
  extract_rbind <- function(attr0) {
    map_reduce(
      dlist0,
      ~.x %>% magrittr::extract2(attr0) %>% as.numeric() %>% t(),
      rbind)
  }
  map_named(names(dlist0[[1]]), extract_rbind)
}


# Collect objects by attributes
collect <- function(l0) {
  list_names <- names(l0[[1]])
  extract_rbind <- function(attr0) {
    map_reduce(l0, ~.x[[attr0]], rbind)
  }
  map_named(list_names, extract_rbind)
}

if_else <- function(test, yes, no) {
  if (test) {return(yes)} else {return(no)}
}

# Extend dim to include vector
dim2 <- function(x) {
  res <- dim(x)
  if_else(is.null(res), length(x), res)
}

# collect_and_reshape <- function(l0) {
#   list_names <- names(l0[[1]])
#   res <- vector("list", length(list_names))
#   for (i in seq_along(list_names)) {
#     res[[i]] <- purrr::map(l0, ~as.numeric(.x[[list_names[i]]])) %>%
#       do.call(rbind, .)
#   }
#   purrr::set_names(res, list_names)
# }
