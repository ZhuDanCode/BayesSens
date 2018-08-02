# Functions to reshape vectors, matrices and arrays
vec_to_mat <- vec_to_arr <- function(v0, dim0) {
  if (missing(dim0)) dim0 <- attr(v0, "dim0")
  array(v0, dim = dim0)
}

mat_to_vec <- arr_to_vec <- function(x) {
  res <- as.numeric(x)
  attr(res, "dim0") <- dim(x)
  res
}

# Functions to convert between list of vectors and vector
listv_to_vec <- function(l0) {
  res <- unlist(l0)
  attr(res, "dim0") <- purrr::map_dbl(l0, length)
  res
}

vec_to_listv <- function(v0, dim0) {
  if (missing(dim0)) dim0 <- attr(v0, "dim0")
  end <- cumsum(dim0)
  start <- head(c(1, end + 1), length(end))
  purrr::map2(start, end, ~v0[.x:.y])
}

# Functions to convert between list of array and vector
slice_list <- function(x) {
  dim0 <- Map(dim2, x)
  res <- unlist(x)
  attr(res, "dim0") <- dim0
  attr(res, "names0") <- names(x)
  res
}

stack_list <- function(x, dim0, names0) {
  if (missing(dim0)) dim0 <- attr(x, "dim0")
  if (missing(names0)) names0 <- attr(x, "names0")

  arr_ind <- dim0 %>% purrr::map_lgl(~length(.x) > 1)
  el_len  <- dim0 %>% purrr::map_dbl(prod)

  res <- vec_to_listv(x, el_len)
  res[arr_ind] <- purrr::map2(
    .x = res[arr_ind],
    .y = dim0[arr_ind],
    .f = ~vec_to_arr(.x, .y)
  )
  res[!arr_ind] %<>% purrr::map(as.vector)
  setNames(res, names0)
}

# Helper functions
#' Extend dim to include vector
dim2 <- function(x) {
  res <- dim(x)
  if_else(is.null(res), length(x), res)
}

if_else <- function(test, yes, no) {
  # Note that 'no' is eager-evaluated
  if (test) {return(yes)} else {return(no)}
}
