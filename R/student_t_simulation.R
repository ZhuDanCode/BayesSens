#' Simulate data from a linear regression model with student-t noise
#' @param n integer; number of datapoints.
#' @param p integer; number of measurements.
#' @param beta (optional) numeric vector; the regression coefficients.
#' It should be of length (p+1) if 'intercept = T', and the intercept
#' should be the first element.
#' @param df positive number; degrees of freedom, Inf is allowed.
#' @param intercept TRUE or FALSE; whether to include an intercept / bias term.
#' @return A list of the data X, y, and the parameters beta, df.
#' @examples
#' data0 <- student_t_data(n = 100, p = 10)
#' data1 <- student_t_data(n = 100, p = 3, beta = c(0.2, -0.3, 0.4), df = 2, intercept = FALSE)
#' @export
student_t_data <- function(n, p, beta, df, intercept = TRUE) {
  num_covariates <- ifelse(intercept, p + 1, p)
  if (missing(beta)) beta <- rnorm(num_covariates)
  if (length(beta) != num_covariates)
    stop(paste("The length of beta does not match the number of covariates.",
               "Did you forget to specify the intercept?\n"))
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  if (intercept) X <- cbind(1, X)

  if (missing(df)) df <- sample(10, 1)
  y <- X %*% beta + rt(n, df = df)
  list(X = X, y = y, beta = beta, df = df)
}
