#'@import Matrix
#'@export
#run = 'all' means run cv.GTV in sequence (i.e. loop through each combination of lam_TV and lam_ridge)
#run = 'stage1' means run cv.GTV for one pair of lam_TV and lam_ridge values. Save the results in a file begining with file_identifier.
#run = 'stage2' read in all the files from stage1, merge them, and finish the cross validation.
#If run = 'all', ignore file_identifier, i, and j.
#If run = 'stage1', then set i to the index of the lam_TV value you want and l to the index of the lam_ridge value. Use file_identifier.
#If run = 'stage2', then ignore i and j.
#Example of file_identifier:
#file_identifier = 'out/E'. Then the results for the ith value of lam_TV and lth value of lam_ridge will be saved in 'out/E_i_l.csv'
# E could mean we used emperical covariance. 
cv.GTV<-function(X,y,Sigma,family='gaussian',nlambda=100,nfolds=5,metric=c("mclr","l2"),run='all',file_identifier='',i=NULL,l=NULL){
  lam_TV = exp(seq(log(0.1),log(100),length.out = 5))
  lam_ridge = sort(exp(seq(log(10e-5),log(1),length.out = 5)),decreasing = T)
  
  metric = match.arg(metric,c("mclr","l2"))
  
  Bt = gen.Bt(Sigma)
  
  if (run == 'stage1') {
    fit_il = Fit_CV(X = X,y = y,Sigma = Sigma,lam_TV = lam_TV[i],lam_ridge = lam_ridge[l],
                    fit=NULL,metric = metric,nfolds=nfolds,family=family,Bt=Bt)
    write.csv(unlist(fit_il),paste(file_identifier,'_',i,'_',l,'_',k,'.csv',sep = ''))
  } else {
    nd = length(lam_TV); nl = length(lam_ridge)
    
    errMat.min = Matrix(NA,nd,nl); errMat.1se = Matrix(NA,nd,nl)
    se.errMat.min = Matrix(NA,nd,nl); se.errMat.min = Matrix(NA,nd,nl)
    lam_1.min = Matrix(NA,nd,nl);lam_1.1se = Matrix(NA,nd,nl)
    
    for (i in 1:length(lam_TV)) {
      for (l in 1:length(lam_ridge)) {
        if (run == 'all') {
          fit_il = Fit_CV(X = X,y = y,Sigma = Sigma,lam_TV = lam_TV[i],lam_ridge = lam_ridge[l],
                          fit=NULL,metric = metric,nfolds=nfolds,family=family,Bt=Bt)
          coef = t(unlist(fit_il))
          errMat.min[i,j] = coef[1]
          errMat.1se[i,j] = coef[2]
          se.errMat.min[i,j] = coef[3]
          se.errMat.min[i,j] = coef[4]
          lam_1.min[i,j] = coef[5]
          lam_1.1se[i,j] = coef[6]
        } else {
          res = read.csv(paste(file_identifier,'_',i,'_',l,'.csv',sep = ''))
        }
      }
    }
    
    best.ind.min = which(errMat.min==min(errMat.min),arr.ind = T)[1,] # best sparse model
    best.ind.1se = which(errMat.1se==min(errMat.1se),arr.ind = T)[1,] # best sparse model
    
    v.min = c(lam_TV[best.ind.min[1]],lam_ridge[best.ind.min[2]],lam_1.min[best.ind.min[1],best.ind.min[2]])
    v.1se = c(lam_TV[best.ind.1se[1]],lam_ridge[best.ind.1se[2]],lam_1.1se[best.ind.1se[1],best.ind.1se[2]])
    names(v.min) = c("lam_TV","lam_ridge","lam_1")
    names(v.1se) = names(v.min)
    
    # Fit a full model
    fit1 = GTV_v2(X = X,y = y,Sigma = Sigma,lam_TV = v.min[1],lam_ridge = v.min[2],lam_1 = v.min[3],family = family,Bt=Bt)
    fit2 = GTV_v2(X = X,y = y,Sigma = Sigma,lam_TV = v.1se[1],lam_ridge = v.1se[2],lam_1 = v.1se[3],family = family,Bt=Bt)
    
    rownames(errMat.min) = paste("lam_TV",1:length(lam_TV),sep = "")
    colnames(errMat.min) = paste("lam_ridge",1:length(lam_ridge),sep = "")
    dimnames(errMat.1se) = dimnames(errMat.min)
    dimnames(se.errMat.min) = dimnames(errMat.min)
    dimnames(se.errMat.1se) = dimnames(errMat.min)
    dimnames(lam_1.min) = dimnames(errMat.min)
    dimnames(lam_1.1se) = dimnames(errMat.min)
    
    lambdas = list(lam_TV=lam_TV,lam_ridge = lam_ridge, lam_1.min = lam_1.min, lam_1.1se = lam_1.1se)
    
    betaMat = cbind(fit1$beta,fit2$beta); colnames(betaMat) = c("cv.min","cv.1se")
    betaMat_std = cbind(fit1$beta_std,fit2$beta_std); colnames(betaMat_std) = colnames(betaMat)
    
    return(list(beta=betaMat,beta_std=betaMat_std,
                optLambda = list(opt.min=v.min,opt.1se=v.1se),metric=metric,
                cverror=list(err.min = errMat.min, err.1se = errMat.1se,
                             se.err.min = se.errMat.min, se.err.1se = se.errMat.1se),
                lambdas = lambdas))
  }
}