summary.tvgarch <- function (object, ...)
{
  loglik <- as.matrix(object$logLik)
  rownames(loglik) <- "Log-likelihood:"
  colnames(loglik) <- ""
  if (!is.null(object$order.g) && object$order.g[1] != 0) {
    cat("\n")
    model <- NULL
    if (length(object$order.g) == 1) {
      model <- paste(model, "TV(", object$order.g, ")-", sep = "") 
    }
    else {
      for (i in 1:length(object$order.g)){
        if (i == 1) model <- paste(model, "TV(", object$order.g[i], sep = "") 
        else if (i == length(object$order.g)) {
          model <- paste(model, ",", object$order.g[i], ")-", sep = "") 
        }
        else model <- paste(model, ",", object$order.g[i], sep = "")
      }
    }
    if (object$order.h[3] != 0) {
      model <- paste(model, "GJR(", object$order.h[3], ")-", sep = "") 
    }
    if (object$order.h[1] != 0) model <- paste(model, "G", sep = "") 
    model <- paste(model, "ARCH(", sep = "") 
    if (object$order.h[1] != 0) model <- paste(model, object$order.h[1], ",", sep = "") 
    model <- paste(model, object$order.h[2], ")", sep = "") 
    if (!is.null(object$xreg)) model <- paste(model, "-X", sep = "") 
    cat("Model:", model, "\n")
    cat("Transition variable:", colnames(object$xtv), "\n")
    cat("\n")
    if (object$turbo == TRUE) {
      object$se.g <- sqrt(diag(vcov.tvgarch(object = object, spec = "tv")))
    }
    estimates.g <- as.matrix(rbind(object$par.g, object$se.g))
    rownames(estimates.g) <- c("Estimate:", "Std. Error:")
    colnames(estimates.g) <- object$names.g 
    if (object$turbo == TRUE) {
      object$se.h <- sqrt(diag(vcov.tvgarch(object = object, spec = "garch")))
    }
    estimates.h <- as.matrix(rbind(object$par.h, object$se.h))
    rownames(estimates.h) <- c("Estimate:", "Std. Error:")
    colnames(estimates.h) <- object$names.h
    estimates <- cbind(estimates.g, estimates.h)
    print(round(estimates, 7))
    print(round(loglik, digits = 4))
    cat("\n")
  }
  if (is.null(object$order.g) || object$order.g[1] == 0) {
    cat("\n")
    model <- NULL
    if (object$order.h[3] != 0) {
      model <- paste(model, "GJR(", object$order.h[3], ")-", sep = "") 
    }
    if (object$order.h[1] != 0) model <- paste(model, "G", sep = "") 
    model <- paste(model, "ARCH(", sep = "") 
    if (object$order.h[1] != 0) model <- paste(model, object$order.h[1], ",", sep = "") 
    model <- paste(model, object$order.h[2], ")", sep = "") 
    if (!is.null(object$xreg)) model <- paste(model, "-X", sep = "") 
    cat("Model:", model, "\n")
    cat("\n")
    if (object$turbo == TRUE) object$se.h <- sqrt(diag(vcov.tvgarch(object = 
                                                                      object)))
    estimates.h <- as.matrix(rbind(object$par.h, object$se.h))
    rownames(estimates.h) <- c("Estimate:", "Std. Error:")
    colnames(estimates.h) <- object$names.h
    print(round(estimates.h, 7))
    print(round(loglik, digits = 4))
    cat("\n")
    
  }
}