vcov.tvgarch <- function(object, spec = NULL, ...)
{
  if (!is.null(object$order.g)){  
    npar.g <- length(object$par.g)
    if (object$turbo == TRUE) {
      if (is.null(spec) || spec == "tv") {
        s <- length(object$order.g)
        object$vcov.g <- matrix(NA, npar.g, npar.g)
        jac.g <- jacobian(func = tvObj, x = object$par.g[-c(1,(1+s+1):(1+2*s))], fixed.par.g = object$par.g[c(1,(1+s+1):(1+2*s))], xtv = object$xtv, opt = object$opt, order.g = object$order.g, 
                          fixed.h = object$h, y = object$y, iter0 = FALSE, flag = 0)
        J.g <- crossprod(jac.g)  
        H.g <-  optimHess(par = object$par.g[-c(1,(1+s+1):(1+2*s))], fn = tvObj, fixed.par.g = object$par.g[c(1,(1+s+1):(1+2*s))], xtv = object$xtv, opt = object$opt, order.g = object$order.g, 
                          fixed.h = object$h, y = object$y, iter0 = FALSE, flag = 1)
        vcov.g <- solve(-H.g) %*% J.g %*% solve(-H.g)
        object$vcov.g[c(2:(s+1),(2*s+2):npar.g),c(2:(s+1),(2*s+2):npar.g)] <- vcov.g
        rownames(object$vcov.g) <- object$names.g
        colnames(object$vcov.g) <- object$names.g
      } 
      if (is.null(spec) || spec == "garch") {
        fit.h <- garchx(y = object$y/sqrt(object$g), order = object$order.h, xreg = object$xreg)
        object$vcov.h <- vcov.garchx(object = fit.h)
      }
    }
    if (is.null(spec)) {
        npar <- length(object$par.h) + npar.g
        object$vcov <- matrix(NA, npar, npar)
        object$vcov[1:npar.g,1:npar.g] <- object$vcov.g
        object$vcov[(npar.g+1):npar,(npar.g+1):npar] <- object$vcov.h
        colnames(object$vcov) <- c(object$names.g, object$names.h)
        rownames(object$vcov) <- c(object$names.g, object$names.h)
        return(object$vcov)
      }
    if (spec == "tv") return(object$vcov.g)
    if (spec == "garch") return(object$vcov.h)
  }
  if (is.null(object$order.g)) return(vcov.garchx(object = object))
}
