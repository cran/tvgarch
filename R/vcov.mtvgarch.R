vcov.mtvgarch <- function (object, spec = c("sigma2", "tv", "garch", "cc"), ...)
{
  spec <- match.arg(spec)
  if (spec != "cc") {
    results <- list()
    m <- ncol(object$y)
    for (i in 1:m) {
      name <- object$names.y[i]
      object.i <- object$Objs[[paste("obj", i, sep = "")]]
      if (is.null(object$order.g) || any(object$order.g == 0)) {
        if (spec == "tv" || spec == "garch") {
          stop ("Only GARCH models have been estimated.")
        }
      }
      if (is.null(object$order.g) || object$order.g[i,1] == 0) {
          spec.i <- "sigma2"
      }
      if (!is.null(object$order.g) && object$order.g[i,1] != 0) {
        spec.i <- spec
      }
      results[[paste("vcov", paste(spec.i), name, sep = "_")]] <- 
        vcov.tvgarch(object = object.i, spec = spec.i)
    }    
  }
  if (spec == "cc") {  
    if (object$turbo == TRUE) {
      jac.dcc <- jacobian(func = dccObj, x = object$par.dcc, 
                          z = object$residuals, sigma2 = object$sigma2, 
                          flag = 0)
      J.dcc <- crossprod(jac.dcc)  
      H.dcc <- optimHess(par = object$par.dcc, fn = dccObj, 
                         z = object$residuals, sigma2 = object$sigma2, flag = 1)
      solHdcc <- solve(-H.dcc)
      object$vcov.dcc <- solHdcc %*% J.dcc %*% solHdcc
    }
    results <- object$vcov.dcc
  }
  return(results)
}
