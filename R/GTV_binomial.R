#'@import glmnet
#'@import Matrix
#'@export
#' This is for GTV binomial
GTV_binomial<-function(XN = XN,y = y,Bt = Bt,delta = delta,lam_ridge = lam_ridge,lam_1TV = lam_1TV,
                       nlambda = nlambda,lambda.min.ratio = lambda.min.ratio,
                       eps = eps, inner_eps = 1e-4, init = NULL, Hessian_bound=F, maxit = 1000){
  
  # Input: XN, y, Bt, weights, delta, lam_ridge, lam_1TV (or nlambda and lambda.min.ratio)
  m = nrow(Bt); n = nrow(XN); p = ncol(XN)
  r = create_BigXD(XN = XN,Bt = Bt,delta = delta,lam_ridge = lam_ridge)
  BigXD = r$BigXD; invB = r$invB; remove(r)
  
  if(is.null(init)){param = c(mean(y),rep(0,p))}else{param = init}
  if (is.null(lam_1TV)) {
    resp = y - mean(y) 
    BigY = c(resp,rep(0,m))
    fit0<-lm(BigY~c(rep(1,n),rep(0,m))-1)
    lammax = max(abs(t(BigXD[,-1])%*%fit0$residuals/nrow(BigXD)))
    lam_1TV = exp(seq(log(lammax),log(lambda.min.ratio*lammax),length.out = nlambda))
  }
  
  pf = c(0,rep(1,(m+p))) # no penalty for an intercept
  coefMat <- Matrix(0,nrow = (p+1),ncol = length(lam_1TV))
  converged = F; iters=c(); 
  
  # loss <-Matrix(0,nrow = maxit,ncol = length(lam_1TV))
  
  # First iteration
  eta = matMult(XN,param[-1])+param[1]
  pi = 1/(1+exp(-eta))
  if(Hessian_bound){wei = rep(1/4,nrow(XN))}else{wei = pi*(1-pi)}
  BigW = c(wei,rep(1,m))
  resp = (y - pi)/wei + eta
  BigY = c(resp,rep(0,m))
  
  for(k in 1:length(lam_1TV)){
    # print(k)
    it = 1
    while(!converged&&it<maxit){
      
      glmfit1 = glmnet(x = BigXD,y = BigY,weights = BigW,lambda = lam_1TV[k],intercept = F,
                       standardize = F,penalty.factor = pf,thresh = inner_eps)
      
      param.new = c(glmfit1$beta[1],invB%*%glmfit1$beta[-1])
      diff = max(abs(param.new-param))
      
      if(diff<eps){converged=T}else{
        param = param.new; it=it+1
        eta = eta = matMult(XN,param[-1])+param[1]; pi = 1/(1+exp(-eta))
        if(!Hessian_bound){wei = pi*(1-pi); BigW = c(wei,rep(1,m))}
        resp = (y - pi)/wei + eta
        BigY = c(resp,rep(0,m))}
    }
    coefMat[,k] = param.new
    iters[k] = it
    converged=F
  }
  
  return(list(coef_std = coefMat,delta = delta,lam_ridge = lam_ridge,lam_1TV = lam_1TV,iters= iters))
}
