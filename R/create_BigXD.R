# #' @export
create_BigXD<-function(XN,Bt,lam_TV,lam_S){
  
  m = nrow(Bt);
  n = nrow(XN);
  p = ncol(XN)
  
  if(lam_S>0){
    BigX = cbind(c(rep(1,n),rep(0,m)),rbind(XN,sqrt(2*n*lam_S)*Bt)) #(n+m) by (p+1)
  }else{
    BigX = cbind(rep(1,n), XN) # when lam_S = 0. Just adding an intercept column
  }
  
  
  ## Generalized inverse of tilde_Gamma(lam_TV) = ginv([lam_TV*Gamma, Ip])
  if(lam_TV>0){
    B = rbind(lam_TV*Bt,Matrix::diag(x = 1,nrow = p,ncol = p)); # (m+p) by p
    invB = MASS::ginv(as.matrix(B)) # p by (m+p) 
  }else{
    B = Matrix::diag(x = 1,nrow = p,ncol = p)
    invB = Matrix::diag(x = 1,nrow = p,ncol = p)
  }
  
  
  ##  D =  [1,    0_{m+p}'
  ##        0_p, tilde_Gamma(lam_TV)]
  if(lam_TV > 0){
    D = cbind(c(1,rep(0,p)),rbind(rep(0,m+p),invB))#(p+1) by (m+p+1)
  }else{
    D =  Matrix::diag(x = 1,nrow = p+1,ncol = p+1)
  } 
  
  
  if(Matrix::nnzero(BigX)/length(BigX)<0.1){ # if X sparse
    D = Matrix(D,sparse = T)
    BigXD = BigX%*%D
  }else{
    BigX = as.matrix(BigX) 
    BigXD = BigX%*%D # 
  }
  
  ## Return XD ((n+m) by (m+p+1)), tilde_Gamma(lam_TV) (p by (m+p))
  return(list(BigXD=BigXD, invB=invB)) 
}
