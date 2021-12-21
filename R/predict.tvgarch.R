predict.tvgarch <- function (object, n.ahead = 10, newxtv = NULL, 
                             newxreg = NULL, newindex = NULL, n.sim = 5000, 
                             as.zoo = TRUE, verbose = FALSE, ...)
{
  n <- length(object$y)
  if (!is.null(object$order.g) && object$order.g[1] != 0) {
      coefs.h <- as.numeric(coef.tvgarch(object = object, spec = "garch"))
    }
  if (is.null(object$order.g) || object$order.g[1] == 0) {
      coefs.h <- as.numeric(coef.tvgarch(object = object, spec = "sigma2"))
    }
  interceptCoef <- coefs.h[1]
  archCoef <- coefs.h[(1+1):(1+object$order.h[2])]
  if (object$order.h[1] != 0) {
      garchCoef <- coefs.h[(1+object$order.h[2]+1):(1+sum(object$order.h[1:2]))]
    }
  if (object$order.h[1] == 0) garchCoef <- NULL
  if (object$order.h[3] != 0) {
      asymCoef <- 
        coefs.h[(1+sum(object$order.h[1:2])+1):(1+sum(object$order.h))]
    }
  if (object$order.h[3] == 0) asymCoef <- NULL
  if (!is.null(object$xreg)) xregCoef <- coefs.h[-(1:(1+sum(object$order.h)))]
  if (is.null(object$xreg)) xregCoef <- NULL
  backcast.values <- list()
  if (max(object$order.h) > 0) {
    backcast.values$innovations <- tail(as.numeric(object$residuals), 
                                        n = max(object$order.h))
    backcast.values$z2 <- backcast.values$innovations^2
    backcast.values$Ineg <- as.numeric(backcast.values$innovations < 0)
    backcast.values$sigma2 <- tail(as.numeric(object$h), 
                                   n = max(object$order.h))
    if (!is.null(object$xreg)) {
      backcast.values$xreg <- tail(as.numeric(object$xreg%*%xregCoef), 
                                   n = max(object$order.h))
    }
  }
  xreg <- NULL
  if (!is.null(object$xreg)) {
      if ((NROW(newxreg) != n.ahead)) {
        stop("NROW(newxreg) is unequal to n.ahead")
      }
      xreg <- cbind(newxreg) %*% xregCoef
    }
  etahat <- coredata(object$residuals)
  draws <- runif(n.ahead*n.sim, min = 0.5+1e-09, 
                 max = length(etahat)+0.5-1e-09)
  draws <- round(draws, digits = 0)
  innovations <- matrix(etahat[draws], n.ahead, n.sim)
  predictions.h <- matrix(NA, n.ahead, n.sim)
  for (i in 1:n.sim) {
      predictions.h[,i] <- garchxSim(n = n.ahead, intercept = interceptCoef, 
                                   arch = archCoef, garch = garchCoef, 
                                   asym = asymCoef, xreg = xreg, 
                                   innovations = innovations[,i], 
                                   backcast.values = backcast.values, 
                                   verbose = TRUE, as.zoo = FALSE)[,"sigma2"]
    }
  if (!is.null(object$order.g) && object$order.g[1] != 0) {
      coefs.g <- as.numeric(coef.tvgarch(object = object, spec = "tv"))
      if (!is.null(newxtv)) {
        if ((NROW(newxtv) != n.ahead)) {
          stop("NROW(newxtv) is unequal to n.ahead")
        }
        newxtv <- as.matrix(newxtv)
        predictions.g <- tvObj(par.g = coefs.g, fixed.par.g = NULL, 
                               xtv = newxtv, opt = object$opt, 
                               order.g = object$order.g, fixed.h = NULL, 
                               y = NULL, iter0 = TRUE, flag = 2)
      }
      if (is.null(newxtv)) predictions.g <- rep(coefs.g[1], n.ahead)
      result <- as.matrix(rowMeans(predictions.h)*predictions.g)
      if (verbose == TRUE) {
        result <- cbind(result, as.vector(predictions.g), predictions.h)
        colnames(result) <- c("p_sigma2", "p_g", paste("p_h", 
                                                       1:NCOL(predictions.h),
                                                       sep = "_"))
      }
    }
  if (is.null(object$order.g) || object$order.g[1] == 0) {
      result <- as.matrix(rowMeans(predictions.h))
      if (verbose == TRUE) {
        result <- cbind(result, predictions.h)
        colnames(result) <- c("p_sigma2", paste("p_h", 1:NCOL(predictions.h), 
                                              sep = "_"))
      }
    }
  if (is.null(newindex)) {
      newindex <- 1:n.ahead
    }
  if (as.zoo == TRUE) {
      result <- zoo(result, order.by = newindex)
    }
  return(result)
}