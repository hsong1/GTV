#'@export
Fit_CV<-function(X,y,Sigma,lam_TV,lam_ridge,fit=NULL,
                 metric = c("mclr","l2"),family='Gaussian',nfolds=5,Bt=NULL){
  # Given delta, lam_ridge, select best lam_1TV
  # Return mlcr or l2 error, lambda_1TV
  n = nrow(X)
  if(is.null(fit)){
    fit1=GTV_v2(X = X,y = y,Sigma = Sigma,lam_TV = lam_TV,lam_ridge = lam_ridge,family = family,Bt=Bt)
  }
  
  metric = match.arg(metric,c("mclr","l2"))
  
  res = matrix(NA,nrow = length(fit1$lambdas$lam_1TV),ncol = nfolds)
  folds = sample(cut(1:nrow(X),breaks = nfolds,labels = FALSE))
  
  for(i in 1:nfolds){
    idx <- which(folds==i)
    # Training and test set
    train_X = X[-idx,];train_y = y[-idx]
    test_X = X[idx,];test_y = y[idx]
    
    # Training 
    fitk = GTV_v2(X = train_X,y = train_y,Sigma=Sigma,
                  lam_TV = lam_TV,lam_ridge = lam_ridge,lam_1TV = fit1$lambdas$lam_1TV,family = family,Bt=Bt)
    
    # Testing
    yhat = cbind(rep(1,nrow(test_X)),test_X)%*%fitk$beta
    
    if(metric=="mclr"){
      yhat_cl = 1*(yhat>0.5)
      res[,i] = apply(test_y - yhat_cl,2,function(x){sum(abs(x))})/nrow(test_X)
    }else{
      res[,i] = apply(test_y - yhat,2,function(x){sum(x^2)})/nrow(test_X) 
    }
  } # for loop end
  
  cvm = apply(res,1,mean); cvsd = apply(res,1,sd)/sqrt(nfolds)
  indmin <- min(which(cvm==min(cvm)))
  ind <-  intersect(which(cvm>=cvm[indmin]+cvsd[indmin]),(1:indmin))
  if(length(ind)==0){ind1se <-  indmin
  } else {
    ind1se <-  max(ind)
  }  
  
  err.min = cvm[indmin]; se.err.min=cvsd[indmin]
  lambda_1TV.min = fit1$lambdas$lam_1TV[indmin]
  err.1se = cvm[ind1se]; se.err.1se=cvsd[ind1se]
  lambda_1TV.1se = fit1$lambdas$lam_1TV[ind1se]
  
  return(list(err.min=err.min,err.1se=err.1se,
              se.err.min =se.err.min,se.err.1se = se.err.1se,
              lambda_1TV.min=lambda_1TV.min[[1]],lambda_1TV.1se=lambda_1TV.1se[[1]]))
}