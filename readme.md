## R package 'BayesSense' (alpha release)

### Package Description
The `BayesSense` performs a comprehensive local sensitivity analysis of MCMC 
output with respect to *all* input parameters (prior hyper-parameters and 
starting values) to assess prior robustness and algorithm convergence in 
models estimated via Gibbs samplers. It implements the forward-mode
Automatic Differentiation and exploits vector/matrix calculus results 
to compute large sets of Jacobian matrices, which contain the first order 
derivatives of all MCMC draws with respect to all input parameters. This allows us to compute a range of sensitivity measures for posterior statistics such as 
prior robustness of posterior means and prediction intervals, as well as 
algorithm convergence measures bases on starting value sensitivities. 

Currently the package supports

1. Gaussian model with independent Normal-Inverse Gamma priors
2. Student-t model with independent Normal-Inverse Gamma-Gamma priors
3. 2-equation model with independent Normal-Inverse Wishart priors


### Installation instruction
```{r}
# install.packages("devtools")
devtools::install_github("ZhuDanCode/BayesSense", build_vignettes = TRUE)
```
