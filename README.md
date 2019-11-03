README
================

### Regression with Graph-based Total-Variation Regularization

### Description

Implementation of regression with graph-based regularization. See Li et
al. (2019) \<arXiv:1803.07658\>

### Installation

Install using **devtools** package:

    devtools::install_github("hsong1/GTV")

### Example

Linear regression with simulated data at particular lambdas,
e.g. \(\lambda_{TV} = 0.1, \lambda_S = 0.1, \lambda_1 = 0.1\).

    library(GTV)
    data(exampleGTV)
    fit0 = with(exampleGTV, gtv(X = X,y = y,Sigma = Sigma,lam_TV = 0.1,lam_S = 0.1,lam_1 = 0.1))
    head(coef(fit0))

Linear regression with simulated data where we choose lambdas via 5-fold
cross-validation.

    set.seed(1234)
    fit1 = with(exampleGTV, cv.gtv(X = X,y = y,Sigma = Sigma,parallel = T))
    head(coef(fit1))
