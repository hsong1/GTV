#'@import Matrix
#'@export
gen.Bt<-function(Sigma){
  # Bt ; m by p matrix where m = # of edges
  p = ncol(Sigma)
  nzro_arr = which(Sigma!=0,arr.ind = T)
  nzrow = nzro_arr[,1]; nzcol = nzro_arr[,2];
  nbrow = length(nzrow)
  Bt = Matrix::Matrix(data = 0,nrow=nbrow,ncol = p)
  tmp = rep(FALSE,nbrow) # track which row is non-empty
  for(i in 1:nbrow){
    if(nzrow[i]<nzcol[i]){
      tmp[i]=TRUE
      # |Sigma_kl|^{1/2} where k = nzrow[i], l = nzrow[i]
      sigma_kl = sqrt(abs(Sigma[nzrow[i],nzcol[i]]))
      Bt[i,nzrow[i]] = sigma_kl; 
      Bt[i,nzcol[i]] = -sign(Sigma[nzrow[i],nzcol[i]])*sigma_kl
    }
  }
  Bt[which(Matrix::rowSums(abs(Bt))!=0),,drop=F]
}