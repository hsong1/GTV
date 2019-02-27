# #' @export
cv.gtv.l1 <- function(X,
                      y,
                      Sigma,
                      lam_TV,
                      lam_S,
                      fit = NULL,
                      metric = "l2",
                      family = 'gaussian',
                      nlambda_1 = 30,
                      nfolds = 5,
                      Bt = NULL) {
  
  # Given lam_TV, lam_S, select best lam_1
  # Return l2 error and  lambda_1
  
  if (is.null(fit)) {
    fit1 = gtv(
      X = X,
      y = y,
      Sigma = Sigma,
      lam_TV = lam_TV,
      lam_S  = lam_S,
      nlambda = nlambda_1,
      family = family,
      Bt = Bt
    )
  }
  
  # metric = match.arg(metric, c("l2", "mclr"))
  
  # Output (res) = where we will save l2 errors
  res = matrix(NA, nrow = length(fit1$lambdas$lam_1), ncol = nfolds)
  # nfolds - folds Cross-Validation
  folds = sample(cut(1:nrow(X), breaks = nfolds, labels = FALSE))
  
  for (i in 1:nfolds) {
    idx <- which(folds == i)
    
    # Training and test set
    train_X = X[-idx, ]
    train_y = y[-idx]
    test_X = X[idx, ]
    test_y = y[idx]
    
    # Training
    fiti = gtv(
      X = train_X,
      y = train_y,
      Sigma = Sigma,
      lam_TV = lam_TV,
      lam_S = lam_S,
      lam_1 = fit1$lambdas$lam_1,
      family = family,
      Bt = Bt
    )
    
    # Testing
    yhat = cbind(rep(1, nrow(test_X)), test_X) %*% fiti$beta
    res[, i] = apply(test_y - yhat, 2, function(x) {sum(x ^ 2)}) / nrow(test_X)
    
  } # for loop end
  
  cvm = apply(res, 1, mean)
  cvsd = apply(res, 1, sd) / sqrt(nfolds)
  indmin <- min(which(cvm == min(cvm)))
  ind <-  intersect(which(cvm >= cvm[indmin] + cvsd[indmin]), (1:indmin))
  
  if (length(ind) == 0) {
    ind1se <-  indmin
  } else {
    ind1se <-  max(ind)
  }
  
  err.min = cvm[indmin]
  se.err.min = cvsd[indmin]
  lambda_1.min = fit1$lambdas$lam_1[indmin]
  err.1se = cvm[ind1se]
  se.err.1se = cvsd[ind1se]
  lambda_1.1se = fit1$lambdas$lam_1[ind1se]
  
  return(
    list(
      lambda_1 = fit1$lambdas$lam_1,
      cvm = cvm,
      cvsd = cvsd,
      err.min = err.min,
      err.1se = err.1se,
      se.err.min = se.err.min,
      se.err.1se = se.err.1se,
      lambda_1.min = lambda_1.min,
      lambda_1.1se = lambda_1.1se
    )
  )
}
