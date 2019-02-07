#'@ import
#'
cv.GTV<-function(X,y,Sigma,family='Guassian',nlambda=100,nfolds=5,metric=c("mclr","l2"),sig=''){
  delta = c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5)
  lam_ridge = sort(c(1,10,10^2,10^3,10^4,10^5,10^6,10^7,10^8,10^9,10^10)/(10^10),
                   decreasing = T)
  
  nd = length(delta); nl = length(lam_ridge)
  
  errMat.min = Matrix(NA,nd,nl); errMat.1se = Matrix(NA,nd,nl)
  se.errMat.min = Matrix(NA,nd,nl); se.errMat.1se = Matrix(NA,nd,nl)
  lam_1TV.min = Matrix(NA,nd,nl);lam_1TV.1se = Matrix(NA,nd,nl)
  
  metric = match.arg(metric,c("mclr","l2"))
  
  for (i in 1:length(delta)) {
    for (l in 1:length(lam_ridge)) {
      res = matrix(NA,nrow = nlambda,ncol = nfolds)
      for (k in 1:nlambda) {
        r = read.csv(paste('out/',sig,'_',i,'_',l,'_',k,'.csv',sep = ''))
        res[k,] = t(r$x)
      }
      cvm = apply(res,1,mean); cvsd = apply(res,1,sd)/sqrt(nfolds)
      indmin <- min(which(cvm==min(cvm)))
      ind <-  intersect(which(cvm>=cvm[indmin]+cvsd[indmin]),(1:indmin))
      if(length(ind)==0){ind1se <-  indmin
      } else {
        ind1se <-  max(ind)
      }
      
      err.min = cvm[indmin]; se.err.min=cvsd[indmin]
      lambda_1TV.min = lam_1TV[indmin]
      err.1se = cvm[ind1se]; se.err.1se=cvsd[ind1se]
      lambda_1TV.1se = lam_1TV[ind1se]
      
      lambda_1TV.min=lambda_1TV.min[[1]]
      lambda_1TV.1se=lambda_1TV.1se[[1]]
      
      errMat.min[i,l] = err.min;
      errMat.1se[i,l] = err.1se;
      se.errMat.min[i,l] = se.err.min;
      se.errMat.1se[i,l] = se.err.1se;
      lam_1TV.min[i,l] = lambda_1TV.min
      lam_1TV.1se[i,l] = lambda_1TV.1se
    }
  }
  
  best.ind.min = which(errMat.min==min(errMat.min),arr.ind = T)[1,] # best sparse model
  best.ind.1se = which(errMat.1se==min(errMat.1se),arr.ind = T)[1,] # best sparse model
  
  v.min = c(delta[best.ind.min[1]],lam_ridge[best.ind.min[2]],lam_1TV.min[best.ind.min[1],best.ind.min[2]])
  v.1se = c(delta[best.ind.1se[1]],lam_ridge[best.ind.1se[2]],lam_1TV.1se[best.ind.1se[1],best.ind.1se[2]])
  names(v.min) = c("delta","lam_ridge","lam_1TV")
  names(v.1se) = names(v.min)
  
  # Fit a full model
  fit1 = GTV(X = X,y = y,Sigma = Sigma,delta = v.min[1],lam_ridge = v.min[2],lam_1TV = v.min[3],family = family,Bt=Bt,p=p,n=n)
  fit2 = GTV(X = X,y = y,Sigma = Sigma,delta = v.1se[1],lam_ridge = v.1se[2],lam_1TV = v.1se[3],family = family,Bt=Bt,p=p,n=n)
  
  rownames(errMat.min) = paste("delta",1:length(delta),sep = "")
  colnames(errMat.min) = paste("lam_ridge",1:length(lam_ridge),sep = "")
  dimnames(errMat.1se) = dimnames(errMat.min)
  dimnames(se.errMat.min) = dimnames(errMat.min)
  dimnames(se.errMat.1se) = dimnames(errMat.min)
  dimnames(lam_1TV.min) = dimnames(errMat.min)
  dimnames(lam_1TV.1se) = dimnames(errMat.min)
  
  lambdas = list(delta=delta,lam_ridge = lam_ridge, lam_1TV.min = lam_1TV.min, lam_1TV.1se = lam_1TV.1se)
  
  betaMat = cbind(fit1$beta,fit2$beta); colnames(betaMat) = c("cv.min","cv.1se")
  betaMat_std = cbind(fit1$beta_std,fit2$beta_std); colnames(betaMat_std) = colnames(betaMat)
  
  return(list(beta=betaMat,beta_std=betaMat_std,
              optLambda = list(opt.min=v.min,opt.1se=v.1se),metric=metric,
              cverror=list(err.min = errMat.min, err.1se = errMat.1se,
                           se.err.min = se.errMat.min, se.err.1se = se.errMat.1se),
              lambdas = lambdas))
}