vcov.tvgarch <- function (object, spec = c("sigma2", "tv", "garch"), ...)
{
  if (!is.null(object$order.g) && object$order.g[1] != 0) { 
    spec <- match.arg(spec)
    npar.g <- length(object$par.g)
    if (object$turbo == TRUE) {
      if (spec != "garch") {
        s <- length(object$order.g)
        object$vcov.g <- matrix(NA_real_, npar.g, npar.g)
        jac.g <- jacobian(func = tvObj, x = object$par.g[-c(1,(1+s+1):(1+2*s))], 
                          fixed.par.g = object$par.g[c(1,(1+s+1):(1+2*s))], 
                          xtv = object$xtv, opt = object$opt, 
                          order.g = object$order.g, 
                          fixed.h = object$h, y = object$y, iter0 = FALSE, 
                          flag = 0)
        J.g <- crossprod(jac.g)  
        H.g <-  optimHess(par = object$par.g[-c(1,(1+s+1):(1+2*s))], fn = tvObj, 
                          fixed.par.g = object$par.g[c(1,(1+s+1):(1+2*s))], 
                          xtv = object$xtv, opt = object$opt, 
                          order.g = object$order.g, 
                          fixed.h = object$h, y = object$y, iter0 = FALSE, 
                          flag = 1)
        solHG <- solve(-H.g)
        vcov.g <- solHG %*% J.g %*% solHG
        object$vcov.g[c(2:(s+1),(2*s+2):npar.g),c(2:(s+1),(2*s+2):npar.g)] <- vcov.g
        rownames(object$vcov.g) <- object$names.g
        colnames(object$vcov.g) <- object$names.g
      } 
      if (spec != "tv") {
        object$vcov.h <- vcov.garchx(object = object$aux.h)
      }
    }
    if (spec == "tv") {
      return(object$vcov.g)
    }
    else if (spec == "garch") {
      return(object$vcov.h)
    }
    else if (spec == "sigma2") {
      npar <- length(object$par.h) + npar.g
      object$vcov <- matrix(NA_real_, npar, npar)
      object$vcov[1:npar.g,1:npar.g] <- object$vcov.g
      object$vcov[(npar.g+1):npar,(npar.g+1):npar] <- object$vcov.h
      colnames(object$vcov) <- c(object$names.g, object$names.h)
      rownames(object$vcov) <- c(object$names.g, object$names.h)
      return(object$vcov)
    }
  }
  if (is.null(object$order.g) || object$order.g[1] == 0) {
    if (object$turbo == TRUE) {
      object$vcov.h <- vcov.garchx(object = object$aux.h)
    }
    return(object$vcov.h)  
  }
}
