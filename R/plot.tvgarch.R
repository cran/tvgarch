plot.tvgarch <- function (x, spec = c("sigma2", "tv", "garch"), ...) 
{
  if (!is.null(x$order.g) && x$order.g[1] != 0) {
    spec <- match.arg(spec)
    if (spec == "sigma2") {
      par(mfrow = c(3,1))
    }
    if (spec == "sigma2" || spec == "tv") {
      garchEst <- tvgarch(x$y, order.g = NULL, order.h = x$order.h, 
                          xreg = x$xreg)
      plot(x$y.index, sqrt(fitted(garchEst)), 
           type = "l", ylab = "", xlab = "", ...)
      lines(x$y.index, sqrt(x$g), col = "blue", ...)
      title(paste("Series: ", paste(colnames(x$y)), 
                  "\n\n GARCH (black) and TV in TV-GARCH (blue)"))
    }
    if (spec == "sigma2") {
      plot(x$y.index, sqrt(x$sigma2), type = "l", 
           ylab = "", xlab = "", ...)
      title("TV-GARCH")
    }
    if (spec == "sigma2" || spec == "garch") {
      plot(x$y.index, sqrt(x$h), type = "l", ylab = "",
           xlab = "", ...)
      title("GARCH in TV-GARCH")
    }
    par(mfrow=c(1,1))
  }
  if (is.null(x$order.g) || x$order.g[1] == 0) {
    plot(x$y.index, sqrt(x$h), type = "l", ylab = "", 
         xlab = "", ...)
    title(paste("Series: ", paste(colnames(x$y)), "\n\n GARCH"))
  }
}