#'@import Matrix
#'@import glmnet
#'@export
GTV<-function(X,y,Sigma,delta,lam_ridge,lam_1TV=NULL,family='Gaussian',nlambda = 100,
              lambda.min.ratio = ifelse(n<p,0.01,1e-04),standardize=T,Bt=NULL,m=NULL,
              p=NULL,n=NULL,B=NULL,invB=NULL,BigX=NULL,BigXD=NULL){
  if (family == 'Gaussian') {
    if (is.null(Bt)) {Bt = gen.Bt(Sigma)} # m by p
    if (is.null(m)) {m = nrow(Bt)}
    if (is.null(p)) {p = ncol(X)}
    if (is.null(n)) {n = nrow(X)}
    if (is.null(B)) {
      # B = [Bt; \deltaIp] (m+p) by p
      B = rbind(Bt,delta*Matrix::diag(x = 1,nrow = p,ncol = p));
    }
    if (is.null(invB)) {invB = MASS::ginv(as.matrix(B))}
    
    # Standardization
    if (standardize) {
      XN = scale(X)
    } else {
      XN = X
    }
    
    # Augmented Matrices
    BigY = c(y,rep(0,m))
    if (is.null(BigX)) {
      BigIX = c(rep(1,n),rep(0,m))
      BigX = cbind(BigIX,rbind(XN,sqrt(2*n*lam_ridge)*Bt)) #(n+m) by (p+1)
    }
    if (is.null(BigXD)) {
      D = cbind(c(1,rep(0,p)),rbind(rep(0,m+p),invB)) #(p+1) by (m+p+1)
      BigXD = BigX%*%D
      BigXD = Matrix(BigXD,sparse = T)
    }
    
    # Fitting 
    pf = c(0,rep(1,(m+p))) # no penalty for an intercept
    
    fit0<-lm(BigY~BigX[,1]-1)
    
    if (is.null(lam_1TV)) {
      lammax = max(abs(t(BigXD[,-1])%*%fit0$residuals/nrow(BigXD)))
      lam_1TV = exp(seq(log(lammax),log(lambda.min.ratio*lammax),length.out = nlambda))
    }
    
    glmfit1 = glmnet(BigXD,BigY,lambda = lam_1TV,intercept = F,
                     standardize = F,penalty.factor = pf)
    
    if(length(glmfit1$lambda)>1){
      a0 = glmfit1$beta[1,]; theta = glmfit1$beta[-1,]
    }else{
      a0 = glmfit1$beta[1]; theta = glmfit1$beta[-1]
    }
    
    beta = invB%*%theta
    coef_std = rbind(a0,beta)
    
    # Unstandardization
    if (standardize) {
      ## Scale variables
      MX = attr(XN,"scaled:center");SX = attr(XN,"scaled:scale")
      ## Unstandardized coefficients
      beta_ustd = apply(beta,2,function(x){x/SX})
      a0_ustd = a0-apply(beta,2,function(x){sum(x*MX/SX)})
      coef = rbind(a0_ustd,beta_ustd)
    }
    
    if(is.null(colnames(X))){rnames=paste("X",1:p,sep="")}else{rnames=colnames(X)}
    rownames(coef_std) = c("(Intercept)",rnames)
    if (standardize) {
      rownames(coef) = c("(Intercept)",rnames)
    }
    lam_1TV = glmfit1$lambda
    lambdas = list(delta=delta,lam_ridge=lam_ridge,lam_1TV=lam_1TV)
    
    if (standardize) {
      return(list(beta=coef,beta_std=coef_std,lambdas=lambdas))
    }
    return(list(beta=beta,beta_std=coef_std,lambdas=lambdas))
  }
  
  else if (family == 'Binomial') {
    return(GTV_logit(X,y,Sigma,delta,lam_ridge,lam_1TV,lambda.min.ratio,nlambda,Bt=Bt,p=p,n=n))
  }
}