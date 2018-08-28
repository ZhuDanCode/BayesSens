#' Import functions
#' @name imports
#' @importFrom stats runif rnorm rgamma setNames dgamma integrate pgamma rt rchisq rWishart
#' @importFrom graphics hist par
#' @importFrom utils head tail setTxtProgressBar txtProgressBar
#' @importFrom magrittr %>% %<>%
#' @importFrom methods as
#' @importFrom graphics lines plot
utils::globalVariables(c(".", "Var1", "Var2"))
NULL


#' Package 'BayesSense'
#' @name BayesSense
#' @title Sensitivity analysis for Bayesian Inference
#' @author
#'   Liana Jacobi [aut, cph], Dan Zhu [aut, cph], Chun Fung Kwok [aut, cre], Kai-Yang Goh [ctb]
#' @description This package performs a comprehensive local sensitivity analysis of MCMC
#' output with respect to all input parameters (prior hyper-parameters and starting values) to assess prior robustness and algorithm convergence in
#' models estimated via Gibbs samplers.
#' It extends the usual Automatic Differentiation to allow for differentiation of
#' matrix decomposition and random samples, and it exploits matrix and vector calculus results
#' to compute large sets of Jacobian matrices that contain the first order
#' derivatives of all MCMC draws with respect to all input parameters.
#' A range of sensitivity measures for posterior statistics can be computed,
#' such as prior robustness of posterior means and prediction intervals,
#' as well as algorithm convergence measures bases on starting value sensitivities.
NULL


#' @useDynLib BayesSense
#' @importFrom Rcpp sourceCpp
NULL
