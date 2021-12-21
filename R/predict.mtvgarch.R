predict.mtvgarch <- function (object, n.ahead = 10, newxtv = NULL, 
                              newxreg = NULL, newindex = NULL, n.sim = 5000, 
                              as.zoo = TRUE, verbose = FALSE, ...)
{
  m <- ncol(object$y)
  n <- nrow(object$y)
  ncol.i <- 0
  if (verbose == FALSE) {
    v <- m
  }
  if (verbose == TRUE) {
    v <- m*(n.sim + 1)
    if (!is.null(object$order.g)) v <- v + colSums(object$order.g != 0)[1]
  }
  if (!is.null(object$xreg)) {
    newxreg <- as.matrix(newxreg)
    if ((NROW(newxreg) != n.ahead)) {
      stop("NROW(newxreg) is unequal to n.ahead")
    }
  }
  results <- matrix(NA, n.ahead, v)
  colnames(results) <- paste("p", 1:v, sep = "_")
  if (object$spillovers == FALSE) {
    for (i in 1:m) {
      object.i <- object$Objs[[paste("obj", i, sep = "")]]
      name <- object$names.y[i]
      if (!is.null(object$xreg)) {
        if (any(object$order.x[i,] != 0)) {
          newxreg.i <- newxreg[,which(object$order.x[i,] == 1)]
        }
        else newxreg.i <- NULL
      }
      result <- predict.tvgarch(object = object.i, n.ahead = n.ahead,
                                newxtv = newxtv, newxreg = newxreg.i,
                                newindex = newindex, n.sim = n.sim,
                                as.zoo = as.zoo, verbose = verbose)
      results[,(ncol.i+1):(ncol.i+NCOL(result))] <- result
      colnames(results)[(ncol.i+1)] <- paste("p_sigma2", name, sep = "_")
      if (verbose == TRUE) {
          if (is.null(object$order.g) || object$order.g[i,1] == 0) {
            colnames(results)[(ncol.i+2):(ncol.i+NCOL(result))] <- 
              c(paste("p_h", 1:n.sim, name, sep = "_"))
          }
        if (!is.null(object$order.g) && object$order.g[i,1] != 0) {
          colnames(results)[(ncol.i+2):(ncol.i+NCOL(result))] <- 
            c(paste("p_g", name, sep = "_"), paste("p_h", 1:n.sim, name, 
                                                   sep = "_"))
        }
      }
      if (i != m) ncol.i <- ncol.i + NCOL(result)
    }
  }
  if (object$spillovers == TRUE) {
    max.p <- max(object$order.h[,1])
    max.q <- max(object$order.h[,2])
    max.r <- max(object$order.h[,3])
    if (any(c(max.q, max.q, max.q) > 1)) {
      stop("Maximum order.h = (1,1,1) with volatility spillovers when predicting
           variances.")
    }
    max.x <- ncol(object$order.x)
    object$par.h[is.na(object$par.h)] <- 0
    intercept.h <- as.vector(object$par.h[,1])
    ARCH <- matrix(0, m, m)
    diag(ARCH) <- object$par.h[,2]
    par.xreg <- object$par.h[,(2+max.q+max.p+max.r):(1+max.q+max.p+max.r+max.x)]
    ARCH[which(object$order.x == 1)] <- par.xreg[which(object$order.x == 1)]
    GARCH <- matrix(0, m, m)
    ASYM <- matrix(0, m, m)
    if (max.p == 1) {
      diag(GARCH) <- object$par.h[,(2+max.q):(1+max.q+max.p)]
    }
    if (max.r == 1) {
      diag(ASYM) <- object$par.h[,(2+max.q+max.p):(1+max.q+max.p+max.r)]
    }
    backcast.values <- list()
    if (max(object$order.h) > 0) {
      backcast.values$innovations <- tail(object$residuals,
                                          n = max(object$order.h))
      backcast.values$z2 <- backcast.values$innovations^2
      backcast.values$Ineg <- backcast.values$innovations < 0
      backcast.values$h <- tail(object$h, n = max(object$order.h))
    }
    etahat <- object$residuals
    if (verbose == TRUE) predictions.h.sim <- NULL
    predictions.h <- matrix(0, n.ahead, m)
    for (i in 1:n.sim) {
      draws <- runif(n.ahead, min = 0.5+1e-09, max = nrow(etahat)+0.5-1e-09)
      draws <- round(draws, digits = 0)
      innovations <- matrix(etahat[draws,], n.ahead, m)
      Ineg <- innovations < 0   
      y2 <- as.vector(backcast.values$z2[1,])
      ht <- as.vector(backcast.values$h[1:max.p,])
      Ineg.y2 <- as.vector((backcast.values$Ineg*y2)[1:max.r,])    
      h <- matrix(NA_real_, n.ahead, m)
      for (t in 1:n.ahead) {
        ht <- intercept.h + ARCH%*%y2 + GARCH%*%ht + ASYM%*%(Ineg.y2)
        h[t,] <- ht
        y2 <- h[t,]*innovations[t,]^2
        Ineg.y2 <- Ineg[t,]*y2
      }
      predictions.h <- predictions.h + h
      if (verbose == TRUE) {
        colnames(h) <- paste(paste("p_h", paste(object$names.y), sep = "_"), i, 
                             sep = "_")
        predictions.h.sim <- cbind(predictions.h.sim, h)
      }
    }
    predictions.h <- predictions.h/n.sim
    if (!is.null(object$order.g) && any(object$order.g != 0)) {
      coefs.g <- object$par.g
      if (!is.null(newxtv)) {
        if ((NROW(newxtv) != n.ahead)) {
          stop("NROW(newxtv) is unequal to n.ahead")
        }
        predictions.g <- matrix(NA_real_, n.ahead, m)
        newxtv <- as.matrix(newxtv)
        for (i in 1:m) {
          predictions.g[,i] <- tvObj(par.g = coefs.g[i,], fixed.par.g = NULL, 
                                 xtv = newxtv, opt = object$opt, 
                                 order.g = object$order.g[i,], fixed.h = NULL, 
                                 y = NULL, iter0 = TRUE, flag = 2)
        }
      }
      if (is.null(newxtv)) {
        predictions.g <- matrix(coefs.g[,1], n.ahead, m, byrow = TRUE)
      }
      colnames(predictions.g) <- paste("p_g", object$names.y, sep = "_")
      results <- predictions.h*predictions.g
      colnames(results) <- paste("p_sigma2", object$names.y, sep = "_")
      if (verbose == TRUE) {
        results <- cbind(results, predictions.g, predictions.h.sim)
      }
    }
    if (is.null(object$order.g) || any(object$order.g == 0)) {
      results <- predictions.h
      colnames(results) <- c(paste("p_sigma2", object$names.y, sep = "_"))
      if (verbose == TRUE) {
        results <- cbind(results, predictions.h.sim)
      }
    }
  }
  if (is.null(newindex)) {
    newindex <- 1:n.ahead
  }
  if (as.zoo == TRUE) {
    results <- zoo(results, order.by = newindex)
  }
  return(results)
}

