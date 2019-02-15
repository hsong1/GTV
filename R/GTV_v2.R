#'@import Matrix
#'@import glmnet
#'@export
GTV_v2<-function(X,y,Sigma,lam_TV,lam_ridge,lam_1=NULL,family=c("gaussian","binomial"),nlambda = 100,
              lambda.min.ratio = ifelse(n<p,0.01,1e-04),weights=NULL,eps=1e-4,inner_eps = 1e-7,
              init=NULL,Hessian_bound=FALSE,maxit=1000){
  
  family = match.arg(family,c("gaussian","binomial"))
  XN = scale(X) # centered and scaled
  Bt = gen.Bt(Sigma)
  
  delta = 1/lam_TV
  if(is.null(lam_1)){lam_1TV=NULL}else{lam_1TV = lam_1*lam_TV}
  
  if (family == "gaussian") {
    if(is.null(weights)){weights = rep(1,nrow(XN))}
    r = GTV_gaussian(XN = XN,y = y,Bt = Bt,weights = weights,
                     delta = delta,lam_ridge = lam_ridge,lam_1TV = lam_1TV,
                     nlambda = nlambda,lambda.min.ratio = lambda.min.ratio,
                     eps = eps)

    coef_std = r$coef_std; lam_1 = r$lam_1TV/lam_TV; iters=NULL
    
    } else if(family == "binomial"){
    
      r = GTV_binomial(XN = XN,y = y,Bt = Bt,delta = delta,lam_ridge = lam_ridge,lam_1TV = lam_1TV,
                       nlambda = nlambda,lambda.min.ratio = lambda.min.ratio,
                       eps = eps, inner_eps = inner_eps, init = init, 
                       Hessian_bound=Hessian_bound,maxit=maxit)
      
      coef_std = r$coef_std; lam_1 = r$lam_1TV/lam_TV; iters=r$iters
  }
    
  # Unstandardization
  
  ## Scale variables
  MX = attr(XN,"scaled:center");SX = attr(XN,"scaled:scale")
  ## Unstandardized coefficients
  a0 = coef_std[1,]; beta = coef_std[-1,,drop=F]
  
  beta_ustd = apply(beta,2,function(x){x/SX})
  a0_ustd = a0-apply(beta,2,function(x){sum(x*MX/SX)})
  coef = rbind(a0_ustd,beta_ustd)
  
  
  if(is.null(colnames(X))){rnames=paste("X",1:p,sep="")}else{rnames=colnames(X)}
  rownames(coef_std) = c("(Intercept)",rnames)
  rownames(coef) = c("(Intercept)",rnames)
  lambdas = list(lam_TV = lam_TV,lam_ridge=lam_ridge,lam_1=lam_1)
  
  return(list(beta=coef,beta_std=coef_std,lambdas=lambdas,iters=iters))
}
