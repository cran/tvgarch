plot.mtvgarch <- function (x, spec = c("sigma2", "tv", "garch"), ...) 
{
  for (i in 1:ncol(x$y)) {
    xObj <- x$Objs[[paste("obj", i, sep = "")]]
    xObj$y <- matrix(xObj$y, length(xObj$y), 1)
    colnames(xObj$y) <- paste(colnames(x$y)[i])
    plot.tvgarch(x = xObj, spec = spec)
  }
}
