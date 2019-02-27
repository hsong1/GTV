#' Solve regression problem with graph-based regularization

#'@import Matrix
#'@import glmnet
#'@useDynLib GTV
#'@export
gtv <- function(X,
                y,
                Sigma,
                lam_TV,
                lam_S,
                lam_1 = NULL,
                family = "gaussian",
                nlambda = 30,
                lambda.min.ratio = ifelse(nrow(X) < ncol(X), 0.01, 1e-04),
                weights = NULL,
                eps = 1e-4,
                init = NULL,
                maxit = 1000,
                Bt = NULL) {
  XN = scale(X) # centered and scaled
  if (is.null(Bt)) {
    Bt = gen.Bt(Sigma)
  } # if Bt is not provided, create new Bt
  if (is.null(weights)) {
    weights = rep(1, nrow(XN))
  }
  
  r = gtv_gaussian(
    XN = XN,
    y = y,
    Bt = Bt,
    weights = weights,
    lam_TV = lam_TV,
    lam_S = lam_S,
    lam_1 = lam_1,
    nlambda = nlambda,
    lambda.min.ratio = lambda.min.ratio,
    eps = eps
  )
  
  lam_1 = r$lam_1
  coef_std = r$coef_std
  iters = NULL
  
  
  # Unstandardization
  
  ## Retrieve scales
  MX = attr(XN, "scaled:center")
  SX = attr(XN, "scaled:scale")
  ## Unstandardized coefficients
  a0 = coef_std[1,]
  beta = coef_std[-1, , drop = F]
  
  beta_ustd = apply(beta, 2, function(x) {
    x / SX
  })
  a0_ustd = a0 - apply(beta, 2, function(x) {
    sum(x * MX / SX)
  })
  coef = rbind(a0_ustd, beta_ustd)
  
  if (is.null(colnames(X))) {
    rnames = paste("X", 1:ncol(X), sep = "")
  } else{
    rnames = colnames(X)
  }
  rownames(coef_std) = c("(Intercept)", rnames)
  rownames(coef) = c("(Intercept)", rnames)
  lambdas = list(lam_TV = lam_TV,
                 lam_S = lam_S,
                 lam_1 = lam_1)
  
  return(structure(
    list(
      beta = coef,
      beta_std = coef_std,
      lambdas = lambdas,
      iters = iters
    ),
    class = "gtv"
  )
  )
}
