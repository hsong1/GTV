#'@import Matrix
#'@import glmnet
# #' @export
gtv_gaussian <-
  function(XN,
           y,
           Bt,
           weights = rep(1, nrow(XN)),
           lam_TV,
           lam_S,
           lam_1 = NULL,
           nlambda = 100,
           lambda.min.ratio = ifelse(n < p, 0.01, 1e-04),
           eps = eps) {
    
    # solve Weighted Least Square problem with GTV penalty.
    
    m = nrow(Bt)
    n = nrow(XN)
    p = ncol(XN)
    
    r = create_BigXD(
      XN = XN,
      Bt = Bt,
      lam_TV = lam_TV,
      lam_S = lam_S
    )
    # Extract outputs
    BigXD = r$BigXD; invB = r$invB; remove(r); gc()
    
    if(lam_S>0) {
      BigY = c(y, rep(0, m))
      BigW = c(weights, rep(1, m))
    } else{
      BigY = y
      BigW = weights
    }
    # Intercept factor for (alpha, theta) in R^{1+(m+p)}. No penalty for an intercept
    pf = c(0, rep(1, (m + p))) 
    
    
    # Default Sequence for lam_1
    if (is.null(lam_1)) {
      fit0 <- lm(BigY ~ rep(1, length(BigY))  - 1)
      lammax = max(abs(t(BigXD[, -1]) %*% fit0$residuals / nrow(BigXD)))
      lam_1 = exp(seq(log(lammax), log(lambda.min.ratio * lammax), length.out = nlambda))
    }
    
    glmfit1 = glmnet(
      x = BigXD,
      y = BigY,
      weights = BigW,
      lambda = lam_1,
      intercept = F,
      standardize = F,
      penalty.factor = pf,
      thresh = eps
    )
    
    if (length(glmfit1$lambda) > 1) {
      a0 = glmfit1$beta[1, ]
      theta = glmfit1$beta[-1, ]
    } else{
      a0 = glmfit1$beta[1]
      theta = glmfit1$beta[-1]
    }
    
    beta = invB %*% theta
    coef_std = rbind(a0, beta)
    return(list(
      coef_std = coef_std,
      lam_TV = lam_TV,
      lam_S = lam_S,
      lam_1 = lam_1
    ))
  }