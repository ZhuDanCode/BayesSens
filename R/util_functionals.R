map_reduce <- function(.x, .f, join, ...) {
  .x %>% purrr::map(.f, ...) %>% do.call(join, .)
}

map_named <- function(.x, .f, .n, ...) {
  if (missing(.n)) .n <- .x
  .x %>% purrr::map(.f, ...) %>% purrr::set_names(.n)
}


# # Memoisation
# # This function should be used with scalar arguments since the checking past-evaluations
# # takes time for large objects (potentially solvable with hash encoding).
# memoise <- function(f) {
#   table0 <- list()
#   function(...) {
#     char_x <- collapse_args(...)
#     res <- table0[[char_x]]
#     if (is.null(res)) {
#       res <- do.call(f, list(...))
#       table0[[char_x]] <<- res
#     }
#     return(res)
#   }
# }
# collapse_args <- function(...) {
#   list(...) %>% as.character() %>% paste(collapse = ",")
# }
