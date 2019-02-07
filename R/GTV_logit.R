#'@import Matrix
#'@export
GTV_logit<-function(X,y,Sigma,delta,lam_ridge,lam_1TV,lambda.min.ratio,nlambda,thresh = 10^-3,Bt=NULL,p=NULL,n=NULL) {
  if(is.null(n)){n = nrow(X)}
  XN = scale(X)
  X_tilda = cbind(rep(1,n), XN)
  
  if (is.null(Bt)) {Bt = gen.Bt(Sigma)}
  m = nrow(Bt)
  if (is.null(p)) {p = ncol(X)}
  B = rbind(Bt,delta*Matrix::diag(x = 1,nrow = p,ncol = p));
  invB = MASS::ginv(as.matrix(B))
  
  # Augmented Matrices
  BigIX = c(rep(1,n),rep(0,m))
  BigX = cbind(BigIX,rbind(XN,sqrt(n*lam_ridge)*Bt)) #(n+m) by (p+1)
  D = cbind(c(1,rep(0,p)),rbind(rep(0,m+p),invB)) #(p+1) by (m+p+1)
  BigXD = BigX%*%D
  BigXD = Matrix(BigXD,sparse = T)
  
  eta = rep(0,ncol(X_tilda))
  eta_old = eta + rep(10*thresh,ncol(X_tilda))
  
  return(GTV_logit_helper(XN,X_tilda,y,Sigma,delta,lam_ridge,lam_1TV,Bt,m,p,n,B,invB,BigX,BigXD,thresh,eta,eta_old))
}
