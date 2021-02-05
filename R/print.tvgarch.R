print.tvgarch <- function (x, ...)
{
  if (!is.null(x$order.g)) {
    loglik <- as.matrix(x$logLik)
    rownames(loglik) <- "Log-likelihood:"
    colnames(loglik) <- ""
    cat("\n")
    cat("\n")
    cat("Date:", x$date, "\n")
    model <- NULL
    model <- paste(model, "TV(", length(x$order.g), ")-", sep = "") 
    if (x$order.h[3] != 0) model <- paste(model, "GJR(", x$order.h[3], ")-", sep = "") 
    if (x$order.h[1] != 0) model <- paste(model, "G", sep = "") 
    model <- paste(model, "ARCH(", sep = "") 
    if (x$order.h[1] != 0) model <- paste(model, x$order.h[1], ",", sep = "") 
    model <- paste(model, x$order.h[2], ")", sep = "") 
    if (!is.null(x$xreg)) model <- paste(model, "-X", sep = "") 
    cat("Model:", model, "\n")
    cat("Method: maximization by parts (algorithm stopped after", x$iter, "iterations) \n")
    cat("No. of observations:", length(x$sigma2),
        "\n")
    cat("\n")  
    cat("TV specification (long-term component): \n")
    cat("Optimization: linearly-constrained \n")
    cat("Transition variable:", colnames(x$xtv), "\n")
    if (!is.null(x$message.g)) cat("Message (constrOptim):", x$message.g, "\n")
    cat("\n")
    if (x$turbo == TRUE) x$se.g <- sqrt(diag(vcov.tvgarch(object = x, spec = "tv")))
    estimates.g <- as.matrix(rbind(x$par.g, x$se.g))
    rownames(estimates.g) <- c("Estimate:", "Std. Error:")
    colnames(estimates.g) <- x$names.g 
    colnames(estimates.g)[1] <- "intercept"
    print(round(estimates.g, 7))
    cat("\n")
    cat("GARCH specification (short-term component): \n")
    cat("Optimization: box-constrained \n")
    if (!is.null(x$message.h)) cat("Message (nlminb):", x$message.h, "\n")
    cat("\n")
    if (x$turbo == TRUE) x$se.h <- sqrt(diag(vcov.tvgarch(object = x, spec = "garch")))
    estimates.h <- as.matrix(rbind(x$par.h, x$se.h))
    rownames(estimates.h) <- c("Estimate:", "Std. Error:")
    colnames(estimates.h) <- x$names.h
    colnames(estimates.h)[1] <- "intercept"
    print(round(estimates.h, 7))
    print(round(loglik, digits = 4))
    cat("\n")
  }
  if (is.null(x$order.g)) print.garchx(x = x)
}
