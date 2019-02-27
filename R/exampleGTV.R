#' Create Example Data Set
#'
#' a block complete graph G that has K connected components and 
#' each connected component is a complete graph with p/K nodes

"exampleGTV"

# set.seed(1114)
# n = 100;
# s = 10;
# r = 0.8
# p = 100
# K = 10
# a = K/p
# 
# ## Population Covariance Matrix
# Sigma = Matrix(0,p,p)
# 
# B = matrix(a*r,p/K,p/K) +a*(1-r)*diag(x = rep(1,p/K))
# Sigma = kronecker(diag(1,K),B)
# # image(Matrix(Sigma))
# 
# library(mvnfast)
# X = mvtnorm::rmvnorm(n = n,mean = rep(0,p),sigma = Sigma)
# beta0 = c(rep(1,s),rep(0,p-s)) # first s components are active
# y = X%*%beta0 + rnorm(n = n, sd = .01)
# 
# exampleGTV = list(X=X,y=y,Sigma=Sigma)
# # usethis::use_data(exampleGTV)
