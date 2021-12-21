logLik.mtvgarch <- function (object, ...)
{
  results <- list()
  m <- ncol(object$y)
  n <- nrow(object$y)
  npar <- 0
  for(i in 1:m){
    object.i <- object$Objs[[paste("obj", i, sep = "")]]
    name <- object$names.y[i]
    if (!is.null(object$order.g) && object$order.g[i,1] != 0) {
      spec <- "tvgarch"
    }
    if (is.null(object$order.g[i,1]) || object$order.g[i,1] == 0){ 
      spec <- "garch"
    }
    results[[paste("logLik", paste(spec), paste(name), sep = "_")]] <- 
      logLik.tvgarch(object = object.i)
    npar <- npar + attr(results[[paste("logLik", paste(spec), paste(name), 
                                       sep = "_")]], "df")
  }
  if (is.null(object$par.dcc)){
    z <- object$residuals
    sigma2 <- object$sigma2
    R <- matrix(1, m, m)
    R[lower.tri(R, FALSE)] <- R[upper.tri(R, FALSE)] <- object$ccc[1,]
    lndetR <- log(det(R))
    invR <- solve(R)
    lf <- -0.5*n*m*log(2*pi) - 0.5*sum(log(sigma2)) - 0.5*n*lndetR - 
      0.5*sum((z %*% invR)*z)
    class(lf) <- "logLik"
    attr(lf, "df") <- npar
    attr(lf, "nobs") <- n
    results[["ccc"]] <- lf
  } 
  if (!is.null(object$par.dcc)) {  
    lf <- object$logLik.dcc
    npar <- npar + 2
    class(lf) <- "logLik"
    attr(lf, "df") <- npar
    attr(lf, "nobs") <- n
    results[["dcc"]] <- lf
  }
  return(results)
}