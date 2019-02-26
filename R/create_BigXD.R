#'@useDynLib GTV
#'@export
create_BigXD<-function(XN,Bt,delta,lam_ridge){
  
  m = nrow(Bt);n = nrow(XN);p = ncol(XN)
  # input: XN, Bt, lam_ridge
  BigX = cbind(c(rep(1,n),rep(0,m)),rbind(XN,sqrt(2*n*lam_ridge)*Bt)) #(n+m) by (p+1)
  
  # input: Bt, delta
  B = rbind(Bt,delta*Matrix::diag(x = 1,nrow = p,ncol = p));
  invB = MASS::ginv(as.matrix(B))
  D = cbind(c(1,rep(0,p)),rbind(rep(0,m+p),invB)) #(p+1) by (m+p+1)
  
  
  if(nnzero(BigX)/length(BigX)<0.1){
    D = Matrix(D,sparse = T)
    BigXD = matMult_sp(BigX,D)
  }else{
    BigX = as.matrix(BigX)
    BigXD = matMult(BigX,D)
  }
  
  return(list(BigXD=BigXD, invB=invB))
}
