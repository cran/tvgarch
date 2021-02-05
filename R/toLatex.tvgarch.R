toLatex.tvgarch <- function(object, digits = 4, ...)
{
  if (!is.null(object$order.g)) {
    if (object$turbo == TRUE) {
      object$se.g <- sqrt(diag(vcov.tvgarch(object = object, spec = "tv")))
      object$se.h <- sqrt(diag(vcov.tvgarch(object = object, spec = "garch")))
    }
    coef.h <- object$par.h
    if (object$order.h[1] != 0 && is.null(object$order.g)) names.h <- c("", paste("y^{2}_{t-",paste(seq(1:object$order.h[2]), sep = ""),"}", sep = ""))
    if (object$order.h[1] != 0 && !is.null(object$order.g)) names.h <- c("", paste("\\dfrac{y^{2}_{t-",paste(seq(1:object$order.h[2]), sep = ""),"}}{g_{t-", paste(seq(1:object$order.h[2]), sep = ""),"}}", sep = ""))
    if (object$order.h[1] != 0 && !is.null(object$order.g)) names.h <- c(names.h, paste("\\widehat{h}_{t-",paste(seq(1:object$order.h[1]), sep = ""),"}", sep = "")) 
    if (object$order.h[1] != 0 && is.null(object$order.g)) names.h <- c(names.h, paste("\\widehat{\\sigma}^2_{t-",paste(seq(1:object$order.h[1]), sep = ""),"}", sep = "")) 
    if (object$order.h[3] != 0) names.h <- c(names.h, paste("y^{2}_{t-",paste(seq(1:object$order.h[3]), sep = ""),"}\\text{I}(y < 0)", sep = "")) 
    if (!is.null(object$xreg)) names.h <- c(names.h, paste("x_{",paste(1:ncol(object$xreg),sep = ""),",t}", sep = ""))
    coefsNames.h <- names.h
    coef.h <- as.numeric(coef.h)
    stderrs <- as.numeric(object$se.h)
    eqtxt.h <- NULL
    for (i in 1:length(coef.h)) {
      ifpluss <- ifelse(i==1, "", " + ")
      eqtxt.h <- paste(eqtxt.h, ifelse(coef.h[i]<0, " - ", ifpluss), 
                     "\\underset{(", format(round(stderrs[i], digits = digits), nsmall = digits), ")}{", format(round(abs(coef.h[i]), 
                     digits = digits), nsmall = digits),"}", coefsNames.h[i], sep = "")
    }
    txtAddEq1 <- " \\\\[1mm]"
    txtAddEq2 <- " \\\\[1mm]"
    eqtxt.h <- paste0("  \\widehat{h}_t & = ", eqtxt.h, "", txtAddEq2, " \n")
    s <- length(object$order.g)
    coef.g <- object$par.g[1:(s+1)]
    coef.G <- object$par.g[-(1:(s+1))]
    coefsNames.g <- ""
    for (i in 1:s) {
      coefsNames.g <- c(coefsNames.g, paste("\\widehat{G}_{", paste(i, sep = ""), "t}", sep = ""))
    }
    coef.g <- as.numeric(coef.g)
    stderrs.g <- as.numeric(object$se.g[1:(s+1)])
    stderrs.tv <- as.numeric(object$se.g[-(1:(s+1))])
    eqtxt.g <- NULL
    for (i in 1:length(coef.g)) {
      ifpluss <- ifelse(i == 1, "", " + ")
      eqtxt.g <- paste(eqtxt.g, ifelse(coef.g[i] < 0, " - ",ifpluss))
      if (i == 1) eqtxt.g <- paste(eqtxt.g, "\\underset{(-)}{", format(round(abs(coef.g[i]), digits = digits), nsmall = digits), "}", sep = "")
      else {
        eqtxt.g <- paste(eqtxt.g, "\\underset{(", format(round(stderrs.g[i], digits = digits), nsmall = digits), ")}{", format(round(abs(coef.g[i]), digits = digits), nsmall = digits), "}",
                          "\\widehat{G}_{", paste(i-1, sep = ""), "t}(", "\\underset{(-)}{", format(round(abs(coef.G[i-1]), digits = digits), nsmall = digits), "},", sep = "")
        for (j in 1:object$order.g[i-1]) {
          eqtxt.g <- paste(eqtxt.g, "\\underset{(", format(round(stderrs.tv[s+j], digits = digits), nsmall = digits), ")}{", format(round(abs(coef.G[s+j]), digits = digits), nsmall = digits), "}", sep = "")
          if (j < object$order.g[i-1]) eqtxt.g <- paste(eqtxt.g, ",", sep = "")
        }
        eqtxt.g <- paste(eqtxt.g, ")", sep = "")
      }
    }
    txtAddEq <- " \\\\[1mm]"
    eqtxt.g <- paste0("  \\widehat{g}_t & = ", eqtxt.g, "", txtAddEq, " \n")
    goftxt <- NULL
    goftxt <- "&"
    iT <- length(object$sigma2)
    goftxt <- paste(goftxt, " \\text{Log-likelihood: }", format(round(as.numeric(object$logLik), digits = digits), nsmall = digits), "\\qquad T = ", iT, 
                    " \n", sep = "")
    cat("\\begin{align*}\n")
    cat(eqtxt.h)
    cat(eqtxt.g)
    cat(goftxt)
    cat("\\end{align*}\n")
  }
  if (is.null(object$order.g)) toLatex.garchx(object = object, digits = digits)
}