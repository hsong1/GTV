#'@import Matrix
#'@import glmnet
#'@export
GTV_gaussian <- function(XN, y, Bt, weights = rep(1,nrow(XN)), delta, lam_ridge, 
                         lam_1TV = NULL, nlambda =100, lambda.min.ratio = ifelse(n<p,0.01,1e-04), eps = eps){

  # solve Weighted Least Square problem with GTV penalty.
  # Input: XN, y, Bt, weights, delta, lam_ridge, lam_1TV (or nlambda and lambda.min.ratio)
  
  m = nrow(Bt); n = nrow(XN); p = ncol(XN)
  r = create_BigXD(XN = XN,Bt = Bt,delta = delta,lam_ridge = lam_ridge)
  BigXD = r$BigXD; invB = r$invB; remove(r)
  BigY = c(y,rep(0,m))
  pf = c(0,rep(1,(m+p))) # no penalty for an intercept
  BigW = c(weights,rep(1,m))
  
  if (is.null(lam_1TV)) {
    fit0<-lm(BigY~c(rep(1,n),rep(0,m))-1)
    lammax = max(abs(t(BigXD[,-1])%*%fit0$residuals/nrow(BigXD)))
    lam_1TV = exp(seq(log(lammax),log(lambda.min.ratio*lammax),length.out = nlambda))
  }
  
  glmfit1 = glmnet(x = BigXD,y = BigY,weights = BigW,lambda = lam_1TV,intercept = F,
                   standardize = F,penalty.factor = pf,thresh = eps)
  
  if(length(glmfit1$lambda)>1){
    a0 = glmfit1$beta[1,]; theta = glmfit1$beta[-1,]
  }else{
    a0 = glmfit1$beta[1]; theta = glmfit1$beta[-1]
  }
  
  beta = invB%*%theta
  coef_std = rbind(a0,beta)
  return(list(coef_std = coef_std,delta = delta,lam_ridge = lam_ridge,lam_1TV = lam_1TV))
}