summary.mtvgarch <- function (object, ...)
{
  for (i in 1:ncol(object$y)) {
    cat("\n")
    cat(colnames(object$y)[i], ": \n", sep = "")
    summary.tvgarch(object = object$Objs[[paste("obj", i, sep = "")]], ...)
  }
  if (is.null(object$par.dcc)) {
    cat("\n")
    cat("Constant Conditional Correlations: \n")
    cat(" \n")
    z <- object$residuals
    nobs <- nrow(object$y)
    m <- ncol(object$y)
    sigma2 <- object$sigma2
    lndetR <- log(det(object$R))
    invR <- solve(object$R)
    lf <- -0.5*nobs*m*log(2*pi) -0.5*sum(log(sigma2)) -0.5*nobs*lndetR -
      0.5*sum((z %*% invR)*z)
    print(round(object$R, digits = 4))
    cat(" \n")
    loglik <- as.matrix(lf)
    rownames(loglik) <- "Log-likelihood:"
    colnames(loglik) <- ""
    print(round(loglik, digits = 4))
  }
  if (!is.null(object$par.dcc)) {
    cat("\n")
    cat("Dynamic Conditional Correlations \n")
    cat(" \n")
    cat(" \n")
    cat("Q matrix \n")
    cat(" \n")
    print(round(object$R, digits = 4))
    cat("\n")
    cat("Coefficients \n")
    cat(" \n")
    if (object$turbo == TRUE) {
      object$se.dcc <- 
        matrix(sqrt(diag(vcov.mtvgarch(object = object, 
                                       spec = "cc"))), 1, 2)
    }
    estimates.dcc <- as.matrix(rbind(object$par.dcc, object$se.dcc))
    rownames(estimates.dcc) <- c("Estimate:", "Std. Error:")
    colnames(estimates.dcc) <- c("alpha", "beta")
    print(round(estimates.dcc,4))
    loglik <- as.matrix(object$logLik.dcc)
    rownames(loglik) <- "Log-likelihood:"
    colnames(loglik) <- ""
    print(round(loglik, digits = 4))
    cat("\n")
  }
}