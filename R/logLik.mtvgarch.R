logLik.mtvgarch <- function(object, ...)
{
  results <- list()
  m <- ncol(object$y)
  nobs <- nrow(object$y)
  npar <- 0
  for(i in 1:m){
    object.i <-object$Objs[[paste("obj", i, sep = "")]]
    if (!is.null(object$order.g) && object$order.g[i,1] != 0) results[[paste("tvgarch", i, sep = "")]] <- logLik.tvgarch(object = object.i)
    if (is.null(object$order.g[i,1]) || object$order.g[i,1] == 0){ 
      object.i$maxpqrpluss1 <- 1
      results[[paste("garch", i, sep = "")]] <- logLik.garchx(object = object.i)
    }
    npar <- npar + length(coef.tvgarch(object = object.i))
  }
  if (is.null(object$par.dcc)){
    z <- residuals.mtvgarch(object = object)
    sigma2 <- fitted.mtvgarch(object = object)$sigma2
    R <- matrix(1, m, m)
    R[lower.tri(R, FALSE)] <- R[upper.tri(R, FALSE)] <- object$ccc[1,]
    lndetR <- log(det(R))
    invR <- solve(R)
    lf <- -0.5*nobs*m*log(2*pi) - 0.5*sum(log(sigma2)) - 0.5*nobs*lndetR - 0.5*sum((z %*% invR)*z)
    attr(lf, "df") <- npar
    attr(lf, "nobs") <- nobs
    results[["ccc"]] <- lf
  } 
  if (!is.null(object$par.dcc)) {  
    lf <- object$logLik.dcc
    npar <- npar + 2
    attr(lf, "df") <- npar
    attr(lf, "nobs") <- nobs
    results[["dcc"]] <- lf
  }
  return(results)
}