#' @export
#' @method coef gtv
coef.gtv <- function(object, std.scale = F, ...) {
  fit = object
  if (std.scale) {
    fitcoef = fit$beta_std
  } else{
    fitcoef = fit$beta
  }
  return(fitcoef)
}

#' @export
#' @method coef cvgtv
coef.cvgtv <- function(object, std.scale = F, ...) {
  fit1 = object$gtv_f.min
  fit2 = object$gtv_f.1se
  if (std.scale) {
    fitcoef = data.frame(coef.min = fit1$beta_std, coef.1se = fit$beta_std)
    colnames(fitcoef)= c("coef.min","coef.1se")
  } else{
    fitcoef = data.frame(coef.min = fit1$beta, coef.1se =fit2$beta)
    colnames(fitcoef)= c("coef.min","coef.1se")
  }
  
  return(fitcoef)
}
