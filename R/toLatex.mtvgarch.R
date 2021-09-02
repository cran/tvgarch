toLatex.mtvgarch <- function (object, digits = 4, ...)
{
  m <- ncol(object$y)
  names.y <- colnames(object$y)
  names.sigma2 <- colnames(object$sigma2)
  if(!is.null(object$xreg) && object$spillovers == FALSE) names.x <- colnames(object$xreg)
  if(!is.null(object$xreg) && object$spillovers == TRUE) names.x <- names.y
  cat("Variance equations \n\n") 
  cat("\\begin{eqnarray} \n")
  for(i in 1:m){
    object.i <- object$Objs[[paste("obj", i, sep = "")]]
    object.i$date <- object$date
    object.i$turbo <- object$turbo
    object.i$trace <- object$trace
    if (!is.null(object$order.g) && object$order.g[i,1] != 0) coefs.h <- as.numeric(object.i$par.h)
    if (is.null(object$order.g) || object$order.g[i,1] == 0){
      coefs.h <- as.numeric(coef.garchx(object = object.i))
      object.i$order.h <- object.i$order
    }
    if(object.i$order.h[1] != 0 && is.null(object.i$order.g)) names.h <- c("", paste("y^{2}_{",paste(i, sep = ""),",t-",
      paste(seq(1:object.i$order.h[2]), sep = ""),"}", sep = ""))
    if(object.i$order.h[1] != 0 && !is.null(object.i$order.g)) names.h <- c("", paste("\\dfrac{y^{2}_{",paste(i, sep = ""),",t-",
      paste(seq(1:object.i$order.h[2]), sep = ""),"}}{\\widehat{g}_{",paste(i, sep = ""),",t-", 
      paste(seq(1:object.i$order.h[2]), sep = ""),"}}", sep = ""))
    if(object.i$order.h[1] != 0 && !is.null(object.i$order.g)) names.h <- c(names.h, paste("\\widehat{h}_{",
      paste(i, sep = ""),",t-",paste(seq(1:object.i$order.h[1]), sep = ""),"}", sep = "")) 
    if(object.i$order.h[1] != 0 && is.null(object.i$order.g)) names.h <- c(names.h, paste("\\widehat{\\sigma}^2_{",
      paste(i, sep = ""),",t-",paste(seq(1:object.i$order.h[1]), sep = ""),"}", sep = "")) 
    if(object.i$order.h[3] != 0) names.h <- c(names.h, paste("\\dfrac{y^{2}_{",paste(i, sep = ""),",t-",
      paste(seq(1:object.i$order.h[3]), sep = ""),"}}{\\widehat{g}_{",paste(i, sep = ""),",t-", 
      paste(seq(1:object.i$order.h[3]), sep = ""),"}}\\text{I}(",paste("y_{",paste(i, sep = ""),",t-",
                                                                       paste(seq(1:object.i$order.h[3]), sep = ""),"}", sep = "")," < 0)", sep = "")) 
    if(!is.null(object$xreg) && object$spillovers == FALSE) names.h <- c(names.h, paste("x_{",paste(1:ncol(object$xreg),sep = ""),"t}", sep = ""))
    if(!is.null(object$xreg) && object$spillovers == TRUE){
      for(x in which(object$order.x[i,]==1)){
        if(object$order.g[x,1] == 0) names.h <- c(names.h, paste("y^{2}_{",paste(x,sep = ""),",t-1}", sep = ""))
        if(object$order.g[x,1] == 1) names.h <- c(names.h, paste("\\dfrac{y^{2}_{",paste(x, sep = ""),",t-1}}{\\widehat{g}_{",paste(x, sep = ""),",t-1}}", sep = ""))
      }
    }  
    coefsNames.h <- names.h
    coefs.h <- as.numeric(coefs.h)
    if(object.i$turbo == TRUE) {
      if (is.null(object.i$order.g) || object.i$order.g[1] == 0) object.i$se.h <- sqrt(diag(vcov.tvgarch(object = object.i)))
      if (!is.null(object.i$order.g) || object.i$order.g[1] != 0) {
        object.i$se.h <- sqrt(diag(vcov.tvgarch(object = object.i, spec = "garch")))
        object.i$se.g <- sqrt(diag(vcov.tvgarch(object = object.i, spec = "tv")))
        s <- length(object.i$order.g)
      }
    }
    stderrs <- as.numeric(object.i$se.h)
    eqtxt.h <- NULL
    for (j in 1:length(coefs.h)) {
      ifpluss <- ifelse(j==1, "", " + ")
      eqtxt.h <- paste(eqtxt.h, ifelse(coefs.h[j]<0, " - ",ifpluss), 
                       "\\underset{(", format(round(stderrs[j], digits=digits), nsmall=digits),")}{",format(round(abs(coefs.h[j]), digits=digits), nsmall=digits),"}",coefsNames.h[j], sep="")
    }
    txtAddEq1 <- " \\\\[1mm]"
    txtAddEq2 <- " \\nonumber \\\\[1mm]"
    if (is.null(object$order.g) || object$order.g[i,1] == 0) eqtxt.h <- paste0("  \\sigma^2_{",paste(i,sep = ""),"t} &=& ", eqtxt.h, "", txtAddEq2, " \n")
    if (!is.null(object$order.g) || object$order.g[i,1] != 0) {
      eqtxt.h <- paste0("  \\widehat{h}_{",paste(i,sep = ""),"t} &=& ", eqtxt.h, "", txtAddEq2, " \n")
      if(!is.null(object.i$order.g)){
        s <- length(object.i$order.g)
        coefs.g <- object.i$par.g[1:(s+1)]
        coefs.tv <- object.i$par.g[-(1:(s+1))]
        coefsNames.g <- ""
        for(j in 1:s){
          coefsNames.g <- c(coefsNames.g, paste("\\widehat{G}_{",paste(i, sep = ""),paste(j, sep = ""),"}", sep = ""))
        }
        coefs.g <- as.numeric(coefs.g)
        stderrs.g <- as.numeric(object.i$se.g[1:(s+1)])
        stderrs.tv <- as.numeric(object.i$se.g[-(1:(s+1))])
        eqtxt.g <- NULL
        for (j in 1:length(coefs.g)){
          ifpluss <- ifelse(j == 1, "", " + ")
          eqtxt.g <- paste(eqtxt.g, ifelse(coefs.g[j] < 0, " - ", ifpluss))
          if (j == 1) eqtxt.g <- paste(eqtxt.g, "\\underset{(-)}{", format(round(abs(coefs.g[j]), digits = digits), nsmall = digits),"}", sep = "")
          if (j > 1) {
            eqtxt.g <- paste(eqtxt.g, "\\underset{(", format(round(stderrs.g[j], digits = digits), nsmall = digits),")}{", format(round(abs(coefs.g[j]), 
                       digits = digits), nsmall = digits),"}","\\widehat{G}_{", paste(j-1, sep = ""),"}(",
                       "\\underset{(-)}{", format(round(abs(coefs.tv[j-1]), digits = digits), nsmall = digits),"};", sep = "")
            for(k in 1:object.i$order.g[j-1]){
              eqtxt.g <- paste(eqtxt.g, "\\underset{(", format(round(stderrs.tv[s+k], digits = digits), nsmall = digits),")}{", 
                         format(round(abs(coefs.tv[s+k]), digits = digits), nsmall = digits), "}", sep = "")
              if(k < object.i$order.g[j-1]) eqtxt.g <- paste(eqtxt.g, ",", sep = "")
            }
            if (colnames(object.i$xtv) == "time") eqtxt.g <- paste(eqtxt.g, "; t/n)", sep = "")
            else eqtxt.g <- paste(eqtxt.g, "; s_{t})", sep = "")
          }
        }
        txtAddEq <- " \\\\[1mm]"
        eqtxt.g <- paste0("  \\widehat{g}_{", paste(i, sep = ""),"t} &=& ", eqtxt.g, "", txtAddEq, " \n")
      }
      else eqtxt.g <- paste0("  \\widehat{g}_{", paste(i, sep = ""),"t} &=& 1 \\\\[1mm] \n")
    }
    goftxt <- NULL
    goftxt <- "&&"
    iT <- length(object.i$sigma2)
    goftxt <- paste(goftxt, " \\text{Log-likelihood: }", format(round(as.numeric(logLik.tvgarch(object = object.i)), digits=digits), nsmall=digits), "\\qquad T = ", iT, " \\nonumber \n", sep = "")
    if( i != m) goftxt <- paste(goftxt, " \\\\ [1mm] \n", sep = "")
    cat(eqtxt.h)
    if (!is.null(object$order.g) && object$order.g[i,1] != 0) cat(eqtxt.g)
    cat(goftxt)
  }
  cat("\\end{eqnarray} \n")
  cat(" \n")
  if (is.null(object$par.dcc)) {
    cat("Constant conditional correlations: \n")
    cat(" \n")
    cat("\\begin{table}[!htbp] \\centering \n")
    cat("\\caption{Constant conditional correlations.} \n") 
  }
  if (!is.null(object$par.dcc)) {
    cat("Dynamic conditional correlations: \n")
    cat(" \n")
    cat("\\begin{table}[!htbp] \\centering \n")
    cat("\\caption{Dynamic conditional correlations (unconditional matrix)} \n") 
  }
  cat("\\label{} \n") 
  cat("\\begin{tabular}{@{\\extracolsep{5pt}}", paste(rep("c", m + 1), sep = ""), "} \n")
  cat("\\\\[-1.8ex] \\hline \n")
  cat("\\hline \\\\[-1.8ex] \n") 
  namesR <- paste(" & ") 
  mtxR <- NULL 
  for (i in 1:m){
    if (i < m) namesR <- paste(namesR, colnames(object$R)[i], " & ")
    if (i == m) namesR <- paste(namesR, colnames(object$R)[i], " \\\\ \n ")
    for (j in 1:m){
      if (j == 1) mtxR <- paste(mtxR, colnames(round(object$R, digits))[i], " & ")
      if (j < m) mtxR <- paste(mtxR, round(object$R, digits)[i,j], " & ")
      if (j == m) mtxR <- paste(mtxR, round(object$R, digits)[i,j], " \\\\ \n ")
    }
  }
  cat(namesR)
  cat(mtxR)
  cat("\\hline \\\\[-1.8ex] \n") 
  cat("\\end{tabular} \n") 
  cat("\\end{table} \n")   
  if (!is.null(object$par.dcc)) {
    cat("Dynamic conditional correlations: \n")
    coefs <- as.numeric(object$par.dcc)
    stderrs <- as.numeric(object$se.dcc)
    cat("\\begin{table}[!htbp] \\centering \n")
    cat("\\caption{} \n") 
    cat("\\label{} \n") 
    cat("\\begin{tabular}{@{\\extracolsep{5pt}} cc} \n")
    cat("\\\\[-1.8ex] \\hline \n")
    cat("\\hline \\\\[-1.8ex] \n") 
    cat(" $ \\alpha $ & $ \\beta $ \\\\ \n")
    cat("\\hline \\\\[-1.8ex] \n")
    cat(paste(round(coefs[1], digits), sep = "")," & ", paste(round(coefs[2], digits), sep = ""),"\\\\ \n")
    cat("(",paste(round(stderrs[1], digits), sep = ""),") & (", paste(round(stderrs[2], digits), sep = ""),")\\\\ \n")
    cat("\\hline \\\\[-1.8ex] \n") 
    cat("\\end{tabular} \n") 
    cat("\\end{table} \n") 
  }
}

