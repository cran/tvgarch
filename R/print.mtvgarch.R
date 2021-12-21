print.mtvgarch <- function (x, ...)
{
  cat("\n")
  cat("**************************** \n")
  cat("* The Volatility Component * \n")
  cat("**************************** \n")
  for (i in 1:ncol(x$y)) {
    cat("\n")
    cat("*** Equation ", i, " (",colnames(x$y)[i], ") ***\n", sep = "")
    print.tvgarch(x = x$Objs[[paste("obj", i, sep = "")]])
  }
  cat("\n")
  cat("***************************** \n")
  cat("* The Correlation Component * \n")
  cat("***************************** \n")
  if (is.null(x$par.dcc)) {
    cat("\n")
    cat("Constant Conditional Correlations: \n")
    cat(" \n")
    z <- x$residuals
    nobs <- nrow(x$y)
    m <- ncol(x$y)
    sigma2 <- x$sigma2
    lndetR <- log(det(x$R))
    invR <- solve(x$R)
    lf <- -0.5*nobs*m*log(2*pi) - 0.5*sum(log(sigma2)) - 0.5*nobs*lndetR - 
      0.5*sum((z %*% invR)*z)
    print(round(x$R, digits = 4))
    cat(" \n")
    loglik <- as.matrix(lf)
    rownames(loglik) <- "Log-likelihood:"
    colnames(loglik) <- ""
    print(round(loglik, digits = 4))
  }
  if (!is.null(x$par.dcc)) {
    cat("\n")
    cat("Dynamic Conditional Correlations \n")
    cat(" \n")
    cat(" \n")
    cat("Q matrix \n")
    cat(" \n")
    print(round(x$R, digits = 4))
    cat("\n")
    cat("Coefficients \n")
    cat(" \n")
    if (x$turbo == TRUE) {
      x$se.dcc <- 
        matrix(sqrt(diag(vcov.mtvgarch(object = x, 
                                       spec = "cc"))), 1, 2)
    }
    estimates.dcc <- as.matrix(rbind(x$par.dcc, x$se.dcc))
    rownames(estimates.dcc) <- c("Estimate:", "Std. Error:")
    colnames(estimates.dcc) <- c("alpha", "beta")
    print(round(estimates.dcc,4))
    loglik <- as.matrix(x$logLik.dcc)
    rownames(loglik) <- "Log-likelihood:"
    colnames(loglik) <- ""
    print(round(loglik, digits = 4))
    cat("\n")
  }
}
  
