## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(autoBayes)
n <- 100
p <- 3
data0 <- gaussian_data(n, p, intercept = T)
str(data0)

## ------------------------------------------------------------------------
res <- gaussian_Gibbs(data0$X, data0$y,
  b_0 = rnorm(p+1), B_0 = pdmatrix(p+1)$Sigma, # add one for intercept
  alpha_0 = 13, delta_0 = 8,
)

## ---- fig.width=7, fig.height=6------------------------------------------
plot_posterior(res)

