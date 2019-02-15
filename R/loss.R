#'@export
loss<-function(param,X,y,Sigma,lam_TV,lam_ridge,lam_1,family=c("gaussian","binomial")){
  family = match.arg(family,c("gaussian","binomial"))
  Bt = gen.Bt(Sigma)
  beta = param[-1]; a0 = param[1]
  
  v= Bt%*%beta
  if(family=="gaussian"){
    l = crossprod(y - X%*%beta-a0)/(2*nrow(X))
  }else{
    eta= X%*%beta+a0
    l = sum(log(1+exp(eta))-y*eta)
  }
  
  Pen = lam_ridge*crossprod(v)+lam_1*lam_TV*sum(abs(v))+lam_1*sum(abs(beta))
  Pen = as.numeric(Pen)
  c(l=l,pen=Pen,loss=l+Pen)
}
