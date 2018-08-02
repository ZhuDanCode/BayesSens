#' Import functions
#' @name imports
#' @importFrom stats runif rnorm rgamma setNames dgamma integrate pgamma rt rchisq
#' @importFrom graphics hist par
#' @importFrom utils tail setTxtProgressBar txtProgressBar
#' @importFrom magrittr %>% %<>%
#' @importFrom methods as
#' @importFrom graphics lines plot
utils::globalVariables(c(".", "Var1", "Var2"))
NULL


#' Package 'BayesSense'
#' @name BayesSense
#' @title Sensitivity analysis for Bayesian Inference
#' @author Liana Jacobi <ljacobi@unimelb.edu.au>,
#' Dan Zhu <dan.zhu@monash.edu>,
#' Chun Fung Kwok <kwokcf@unimelb.edu.au>,
#' Kai-Yang Goh <kaiyanggoh@hotmail.com>
#' @description This package computes local derivatives of Bayesian outputs, e.g.
#' MCMC posterior estimates and performance measures, via Automatic differentiation.
NULL


#' @useDynLib BayesSense
#' @importFrom Rcpp sourceCpp
NULL
