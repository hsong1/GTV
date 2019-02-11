#'@export
GTV_logit_helper<-function(XN,X_tilda,y,Sigma,delta,lam_ridge,lam_1TV,Bt,m,p,n,B,invB,BigX,BigXD,thresh,eta,eta_old) { 
  it = 1;
  while (max(abs(eta - eta_old)) > thresh) {
    eta_old = eta
    m_ = (rep(1,nrow(X_tilda))+exp(-(X_tilda %*% eta)))^-1
    mu0 = 4*(y-m_) + (X_tilda %*% eta)
    out = GTV(X=XN,y=mu0,Sigma=Sigma,delta=delta,lam_ridge=4*lam_ridge,lam_1TV=4*lam_1TV,
              standardize=F,Bt=Bt,m=m,p=p,n=n,B=B,invB=invB,BigX=BigX,BigXD=BigXD)
    eta = out$beta_std
    it = it+1
  }
  
  if (ncol(out$beta_std) > 1) {
    a0 = out$beta_std[1,]
  } else {
    a0 = out$beta_std[1]
  }
  
  ## Scale variables
  MX = attr(XN,"scaled:center");SX = attr(XN,"scaled:scale")
  ## Unstandardized coefficients
  beta_ustd = apply(out$beta,2,function(x){x/SX})
  a0_ustd = a0-apply(out$beta,2,function(x){sum(x*MX/SX)})
  out$beta = rbind(a0_ustd,beta_ustd)
  rownames(out$beta) = rownames(out$beta_std)
  out$iter = it
  return(out)
}

