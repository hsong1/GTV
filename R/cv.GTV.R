#' Cross-validation for gtv
#'
#' Does k-fold cross-validation for gtv
#'@import Matrix
#'@useDynLib GTV
#'@export

cv.gtv <- function(X,y,Sigma,
                   family='gaussian',
                   nlambdas = c(5,5,30),
                   nfolds=5,
                   metric="l2",
                   lambdas.seq = NULL,
                   parallel = FALSE,
                   nCores =  detectCores()-1){
  
  nlambda_TV = nlambdas[1]
  nlambda_S = nlambdas[2]
  nlambda_1 = nlambdas[3]
  
  if(is.null(lambdas.seq)){
    lam_TV = sort(exp(seq(log(10e-5),log(1),length.out = nlambda_TV)),decreasing = T)
    lam_S  = sort(exp(seq(log(10e-5),log(1),length.out = nlambda_S)),decreasing = T)
  }else{
    lam_TV = lambdas.seq$lam_TV
    lam_S  = lambdas.seq$lam_S
    nlambda_TV = length(lam_TV)
    nlambda_S = length(lam_S)
  }
  
  # metric = match.arg(metric,c("l2","mclr"))
  # family = match.arg(family,c("gaussian", "binomial"))
  
  Bt = gen.Bt(Sigma)

  # Initialization

  errMat.min = Matrix(NA, nlambda_TV, nlambda_S)
  errMat.1se = Matrix(NA, nlambda_TV, nlambda_S)
  se.errMat.min = Matrix(NA, nlambda_TV, nlambda_S)
  se.errMat.1se = Matrix(NA, nlambda_TV, nlambda_S)
  lam_1.min = Matrix(NA, nlambda_TV, nlambda_S)
  lam_1.1se = Matrix(NA, nlambda_TV, nlambda_S)
  
  lambdas = expand.grid(lam_TV,lam_S)
  
  if(!parallel){
    nCores = 1
  }
  cl = makeCluster(nCores)
  doParallel::registerDoParallel(cl)
  
  res = foreach(l=1:nrow(lambdas),.combine = "cbind")%dopar%{
    
    fit_l = cv.gtv.l1(
      X = X,
      y = y,
      Sigma = Sigma,
      lam_TV = lambdas$Var1[l],
      lam_S = lambdas$Var2[l],
      fit = NULL,
      metric = metric,
      nfolds = nfolds,
      nlambda_1 = nlambda_1,
      family = family,
      Bt = Bt
    )
    with(fit_l, c(err.min, err.1se,se.err.min,se.err.1se,lambda_1.min,lambda_1.1se))
    
  }
  stopCluster(cl)
  
  # Retrieve Outputs
  errMat.min = matrix(res[1,],nrow=length(lam_TV),ncol = length(lam_S))
  errMat.1se = matrix(res[2,],nrow=length(lam_TV),ncol = length(lam_S))
  se.errMat.min = matrix(res[3,],nrow=length(lam_TV),ncol = length(lam_S))
  se.errMat.1se = matrix(res[4,],nrow=length(lam_TV),ncol = length(lam_S))
  lam_1.min = matrix(res[5,],nrow=length(lam_TV),ncol = length(lam_S))
  lam_1.1se = matrix(res[6,],nrow=length(lam_TV),ncol = length(lam_S))
  
  rownames(errMat.min) = lam_TV
  colnames(errMat.min) = lam_S
  dimnames(errMat.1se) = dimnames(errMat.min)
  dimnames(se.errMat.min) = dimnames(errMat.min)
  dimnames(se.errMat.1se) = dimnames(errMat.min)
  dimnames(lam_1.min) = dimnames(errMat.min)
  dimnames(lam_1.1se) = dimnames(errMat.min)
  
  best.ind.min = which(errMat.min == min(errMat.min), arr.ind = T)[1, ] # best sparse model
  best.ind.1se = which(errMat.1se == min(errMat.1se), arr.ind = T)[1, ] # best sparse model
  
  v.min = c(lam_TV[best.ind.min[1]],
            lam_S[best.ind.min[2]],
            lam_1.min[best.ind.min[1], best.ind.min[2]])
  v.1se = c(lam_TV[best.ind.1se[1]], 
            lam_S[best.ind.1se[2]],
            lam_1.1se[best.ind.1se[1], best.ind.1se[2]])
  
  names(v.min) = c("lam_TV", "lam_S", "lam_1")
  names(v.1se) = names(v.min)
  
  # Fit a full model
  fit1 = gtv(
    X = X,
    y = y,
    Sigma = Sigma,
    lam_TV = v.min[1],
    lam_S = v.min[2],
    lam_1 = v.min[3],
    family = family,
    Bt = Bt
  )
  fit2 = gtv(
    X = X,
    y = y,
    Sigma = Sigma,
    lam_TV = v.1se[1],
    lam_S = v.1se[2],
    lam_1 = v.1se[3],
    family = family,
    Bt = Bt
  )
  
  lambdas = list(
    lam_TV = lam_TV,
    lam_S = lam_S,
    lam_1.min = lam_1.min,
    lam_1.1se = lam_1.1se
  )
  
  betaMat = cbind(fit1$beta, fit2$beta)
  colnames(betaMat) = c("cv.min", "cv.1se")
  betaMat_std = cbind(fit1$beta_std, fit2$beta_std)
  colnames(betaMat_std) = colnames(betaMat)
  
  
  return(structure(
    list(
      beta = betaMat,
      beta_std = betaMat_std,
      optLambda = list(opt.min = v.min, opt.1se = v.1se),
      metric = metric,
      cverror = list(
        err.min = errMat.min,
        err.1se = errMat.1se,
        se.err.min = se.errMat.min,
        se.err.1se = se.errMat.1se
      ),
      lambdas = lambdas,
      gtv_f.min = fit1,
      gtv_f.1se = fit2
    ),
    class = "cvgtv"
  ))
}
