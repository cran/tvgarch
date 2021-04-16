predict.mtvgarch <- function (object, n.ahead = 10, newxtv = NULL, newxreg = NULL, newindex = NULL, ...)
{
  results <- NULL
  m <- ncol(object$y)
  n <- nrow(object$y)
  innovations <- matrix(rnorm(n.ahead*m), n.ahead, m)
  for(i in 1:m){
    object.i <- object$Objs[[paste("obj", i, sep = "")]]
    if (!is.null(object$order.g) && object$order.g[i,1] != 0) {
      coefs <- as.numeric(object.i$par.h)
      interceptCoef <- coefs[1]
      archCoef <- coefs[(1+1):(1+object.i$order.h[2])]
      garchCoef <- NULL
      if (object.i$order.h[1] != 0) garchCoef <- coefs[(1+object.i$order.h[2]+1):(1+sum(object.i$order.h[1:2]))]
      asymCoef <- NULL
      if (object.i$order.h[3] != 0) asymCoef <- coefs[(1+sum(object.i$order.h[1:2])+1):(1+sum(object.i$order.h))]
      if (sum(object$order.x[i,]) != 0) xregCoef <- coefs[(1+sum(object.i$order.h)+1):(1+sum(object.i$order.h)+sum(object$order.x[i,]))]
    }
    if (is.null(object$order.g) || object$order.g[i,1] == 0) {
      coefs <- as.numeric(coef.garchx(object = object.i))
      interceptCoef <- coefs[1]
      archCoef <- NULL
      if (object.i$archK > 0) {
        archCoef <- rep(0, object.i$archOrder)
        archCoef[object.i$arch] <- coefs[object.i$archIndx]
      }
      garchCoef <- NULL
      if (object.i$garchK > 0) {
        garchCoef <- rep(0, object.i$garchOrder)
        garchCoef[object.i$garch] <- coefs[object.i$garchIndx]
      }
      asymCoef <- NULL
      if (object.i$asymK > 0) {
        asymCoef <- rep(0, object.i$asymOrder)
        asymCoef[object.i$asym] <- coefs[object.i$asymIndx]
      }
      if (object.i$xregK > 0) {
        xregCoef <- coefs[object.i$xregIndx]
      } 
    }
    backcast.values <- NULL
    xreg <- NULL
    predictions.h <- matrix(NA, n.ahead, 1)
    if (object$spillovers == FALSE) {
      if (!is.null(object.i$xreg)) {
        if ((NROW(newxreg) != n.ahead)) {
          stop("NROW(newxreg) is unequal to n.ahead")
        }
        xreg <- cbind(newxreg) %*% xregCoef
      }
    }
    if (object$spillovers == TRUE) {
      y2.i <- as.matrix(rbind(colMeans(innovations^2), as.matrix(innovations[-n.ahead,]^2)))[,which(object$order.x[i,] == 1)]
      xreg <- cbind(y2.i) %*% xregCoef
    }
    predictions.h <- garchxSim(n.ahead, intercept = interceptCoef, arch = archCoef, garch = garchCoef, asym = asymCoef,
                               xreg = xreg, innovations = innovations[,i], backcast.values = backcast.values, 
                               verbose = TRUE, as.zoo = FALSE)[,"sigma2"]
    if (!is.null(object$order.g) && object$order.g[i,1] != 0) { 
      coefs.g <- as.numeric(coef.tvgarch(object = object.i, spec = "tv"))
      if (!is.null(newxtv)) {
        newxtv <- as.matrix(newxtv)
        predictions.g <- tvObj(par.g = coefs.g, fixed.par.g = NULL, xtv = newxtv, opt = object$opt, order.g = object.i$order.g, fixed.h = NULL, y = NULL, iter0 = TRUE, flag = 2)
      }
      if (is.null(newxtv)) predictions.g <- rep(object$g[n,i], n.ahead)
      result.i <- as.matrix(predictions.h*predictions.g)
    }
    if (is.null(object$order.g) || object$order.g[i,1] == 0) { 
      result.i <- as.matrix(predictions.h)
    }
    colnames(result.i) <- paste("sigma2.", i, sep = "")
    results <- cbind(results, result.i)
  }
  if (is.null(newindex)) {
    newindex <- 1:n.ahead
  }
  results <- zoo(results, order.by = newindex)
  return(results)
}