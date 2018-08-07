#' Check what sensitivities are available
#' @param res Output from *_AD function.
#' @param display T or F; if T, print out summary.
#' @export
available_sensitivity <- function(res, display = T) {
  # Names
  start_with <- purrr::partial(grepl, pattern = "d_")
  var_name <- names(res)
  d_output <- var_name[purrr::map_lgl(var_name, start_with)]
  el <- res[[d_output[1]]]
  d_input <- names(el)

  # Dimension
  output <- setdiff(var_name, d_output)
  output_dim <- purrr::map_dbl(res[output], ncol) %>% setNames(output)
  iter_num <- nrow(res[[output[1]]])
  d_output_dim <- output_dim %>% setNames(d_output)

  var <- output[1]
  d_var <- d_output[purrr::map_lgl(d_output, ~grepl(.x, pattern = var))]
  d_input_dim <- purrr::map_dbl(res[[d_var]], ncol) / output_dim[var]

  # Display
  if (display) {
    cat("Numerator: ", sprintf("%s(%d)", d_output, as.numeric(d_output_dim)), "\n")
    cat("Denominator: ", sprintf("%s(%d)", d_input, as.numeric(d_input_dim)), "\n")
    cat(sprintf("Number of iteration: %d", iter_num), "\n")
  }

  invisible(list(
    numerator = d_output, denominator = d_input,
    output_dim = d_output_dim, input_dim = d_input_dim,
    num_iter = iter_num))
}


#' Retrieve a sensitivitiy from autodiff output
#' @param res Output from *_AD function.
#' @param numerator Character string; the numerator from `available_sensitivity`.
#' @param denominator Character string; the denominator from `available_sensitivity`.
#' @param reshape T or F; if T, reshape the result into an array.
#' @export
get_sensitivity <- function(res, numerator, denominator, reshape = T) {
  el <- res[[numerator]][[denominator]]
  if (!reshape) return(el)

  summary0 <- available_sensitivity(res, display = F)
  array(el, dim = c(summary0$num_iter,
                    summary0$output_dim[numerator],
                    summary0$input_dim[denominator]))
}
