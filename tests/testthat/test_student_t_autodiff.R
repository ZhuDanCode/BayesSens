context("student-t: Check Auto-differentiation results are consistent with Numerical differentiation")
library(BayesSense)

# ===== Helper functions =================================================
# This function is needed to align with the implementation in the Autodiff function.
test_fun <- function(X, y, b_0, B_0, alpha_0, delta_0, nu,
                     init_beta, init_sigma, num_steps = 1e4) {
  if (missing(init_sigma))
    init_sigma <- 1 / sqrt(rgamma(1, alpha_0 / 2, delta_0 / 2))
  if (missing(init_beta))
    init_beta <- MASS::mvrnorm(1, b_0, B_0)

  # Initialisation
  n <- length(y)
  alpha_1 <- alpha_0 + n
  nu_1 <- nu + 1
  beta_g <- init_beta
  sigma_g <- init_sigma
  inv_B_0 <- solve(B_0)
  inv_B_0_times_b_0 <- inv_B_0 %*% b_0
  beta_res <- matrix(0, num_steps, length(b_0))
  sigma_res <- numeric(num_steps)

  pb <- txtProgressBar(1, num_steps, style = 3)
  for (i in 1:num_steps) {
    #Update lambda
    nu_2 <- nu + sigma_g^(-2) * (y - X %*% beta_g) ^ 2
    Lambda <- rgamma(n, nu_1 / 2, nu_2 / 2)

    # Update beta
    XTL <- A_times_diag_v0(t(X), Lambda)
    XTLX <- XTL %*% X
    XTLy <- XTL %*% y
    B_g <- solve(sigma_g^(-2) * XTLX + inv_B_0)
    b_g <- B_g %*% (sigma_g^(-2) * XTLy + inv_B_0_times_b_0)
    # beta_g <- MASS::mvrnorm(1, b_g, B_g)
    beta_g <- b_g + t(chol(B_g)) %*% rnorm(length(b_g))

    # Update sigma
    delta_g <- delta_0 + A_times_diag_v0(t(y - X %*% beta_g), Lambda) %*% (y - X %*% beta_g)
    sigma_g <- 1 / sqrt(rgamma(1, alpha_1 / 2, delta_g / 2))

    # Keep track
    beta_res[i,] <- beta_g
    sigma_res[i] <- sigma_g
    setTxtProgressBar(pb, i)
  }

  # Tidy format
  list(sigma = sigma_res, beta = beta_res)
}

# num_diff_pipe :: list <configuration> -> list <matrices>
num_diff_pipe <- function(config_list, base_case, h) {
  config_list %>%
    purrr::map(~num_diffs(.x, base_case, h)) %>%
    tidy_list() %>%
    purrr::map(t)
}

# num_diffs :: list -> Obj -> list <summary statistics>
num_diffs <- function(new_parameters, base_case, h) {
  set.seed(123)
  new_case <- do.call(test_fun, new_parameters)
  list(d_beta = num_diff(new_case$beta, base_case$beta, h),
       d_sigma2 = num_diff(new_case$sigma^2, base_case$sigma^2, h))
}

# Test if numerical differentiation and auto-differentiation give
# a difference less than a threshold.
test_diff <- function(num_diff, auto_diff) {
  x <- unlist(num_diff)
  y <- unlist(auto_diff)
  print(cbind(num_diff = x, auto_diff = y))
  expect_true(all(abs(x - y) < 1e-5))
}


# ===== Actual testing ===================================================
testthat::test_that("Testing student-t autodiff", {
  # skip_on_cran()
  skip_if_not(Sys.getenv("run_all_tests", TRUE))
  set.seed(123);
  n <- 1000
  p <- 5
  data0 <- student_t_data(n, p, intercept = FALSE)

  set.seed(123);
  b_0 <- rnorm(p)
  B_0 <- pdmatrix(p)$Sigma
  h <- 1e-8

  default_parameters <- list(
    data0$X, data0$y,
    b_0 = b_0, B_0 = B_0, alpha_0 = 3, delta_0 = 5, nu = 20
  )

  set.seed(123)
  base_case <- do.call(test_fun, default_parameters)
  set.seed(123)
  res <- do.call(student_t_AD, append(default_parameters, list(num_steps = 1000)))

  #================== Check sensitivity of beta and sigma wrt b0 ======================
  print("Sensitivity of beta and sigma wrt b0")
  num_diff_sense <- 1:p %>%
    # Generate Configurations
    purrr::map(function(i) {
      new_parameters <- default_parameters
      new_parameters$b_0[i] <- new_parameters$b_0[i] + h
      new_parameters
    }) %>%
    num_diff_pipe(base_case, h)
  auto_diff_sense <- list(
    d_beta = auto_diff(res, "d_beta", "d_b0"),
    d_sigma2 = auto_diff(res, "d_sigma2", "d_b0")
  )
  test_diff(num_diff_sense, auto_diff_sense)

  #================== Check sensitivity of beta and sigma wrt B_0 =====================
  print("Sensitivity of beta and sigma wrt B0")
  num_diff_sense <- expand.grid(1:p, 1:p) %>%
    {purrr::map2(.[,1], .[,2], function(i, j) {
      new_parameters <- default_parameters
      new_parameters$B_0[i, j] <- new_parameters$B_0[i, j] + h
      if (i != j) {
        new_parameters$B_0[j, i] <- new_parameters$B_0[j, i] + h
      }
      new_parameters
    })} %>%
    num_diff_pipe(base_case, h)
  auto_diff_sense <- list(
    beta = auto_diff(res, "d_beta", "d_B0"),
    sigma = auto_diff(res, "d_sigma2", "d_B0")
  )

  test_diff(num_diff_sense, auto_diff_sense)

  #================== Check sensitivity of beta and sigma wrt alpha_0 =================
  print("Sensitivity of beta and sigma wrt alpha_0")
  new_parameters <- default_parameters
  new_parameters$alpha_0 <- new_parameters$alpha_0 + h

  num_diff_sense <- list(new_parameters) %>% num_diff_pipe(base_case, h)

  auto_diff_sense <- list(
    beta = auto_diff(res, "d_beta", "d_alpha0"),
    sigma = auto_diff(res, "d_sigma2", "d_alpha0")
  )

  test_diff(num_diff_sense, auto_diff_sense)

  #================== Check sensitivity of beta and sigma wrt delta_0 =================
  print("Sensitivity of beta and sigma wrt delta_0")
  new_parameters <- default_parameters
  new_parameters$delta_0 <- new_parameters$delta_0 + h

  num_diff_sense <- list(new_parameters) %>% num_diff_pipe(base_case, h)

  auto_diff_sense <- list(
    beta = auto_diff(res, "d_beta", "d_delta0"),
    sigma = auto_diff(res, "d_sigma2", "d_delta0")
  )

  test_diff(num_diff_sense, auto_diff_sense)
})
