#'@export
loss<-function(param,X,y,Sigma,delta,lam_ridge,lam_1TV,family=c("gaussian","binomial")){
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
  
  Pen = as.numeric(lam_ridge*crossprod(v)+lam_1TV*(sum(abs(v))+delta*sum(abs(beta))))
  c(l=l,pen=Pen,loss=l+Pen)
}
