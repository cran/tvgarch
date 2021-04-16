vcov.mtvgarch <- function(object, ...)
{
  results <- list()
  m <- ncol(object$y)
  nobs <- nrow(object$y)
  for(i in 1:m){
    if (!is.null(object$order.g) && object$order.g[i,1] != 0) results[[paste("tvgarch", i, sep = "")]] <- vcov.tvgarch(object = object$Objs[[paste("obj", i, sep = "")]])
    if (is.null(object$order.g[i,1]) || object$order.g[i,1] == 0){ 
      object$Objs[[paste("obj", i, sep = "")]]$maxpqrpluss1 <- 1
      results[[paste("garch", i, sep = "")]] <- vcov.garchx(object = object$Objs[[paste("obj", i, sep = "")]])
    }
  }
  if (!is.null(object$par.dcc)) {  
    if (object$turbo == TRUE) {
      jac.dcc <- jacobian(func = dccObj, x = object$par.dcc, z = object$residuals, sigma2 = object$sigma2, flag = 0)
      J.dcc <- crossprod(jac.dcc)  
      H.dcc <- optimHess(par = object$par.dcc, fn = dccObj, z = object$residuals, sigma2 = object$sigma2, flag = 1)
      object$vcov.dcc <- solve(-H.dcc) %*% J.dcc %*% solve(-H.dcc)
    }
    results[["dcc"]] <- object$vcov.dcc
  }
  return(results)
}
