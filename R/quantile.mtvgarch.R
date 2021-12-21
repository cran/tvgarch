quantile.mtvgarch <- function (x, probs = 0.025, type = 7, as.zoo = TRUE, 
                               ...)
{
  m <- ncol(x$y)
  n <- nrow(x$y)
  ncol.i <- 0
  iCols <- length(probs)
  if (iCols == 1) {
    results <- matrix(NA, n, m)
  }
  if (iCols > 1) {
    results <- matrix(NA, n, m*iCols)
    colnames(results) <- paste("q", 1:ncol(results), sep = "_")
  }  
  for(i in 1:m){
    result.i <- quantile.tvgarch(x = x$Objs[[paste("obj", i, sep = "")]], 
                                 probs = probs, names = TRUE, type = type, 
                                 as.zoo = as.zoo)
    name <- x$names.y[i]
    results[,(ncol.i+1):(ncol.i+iCols)] <- result.i
    if (iCols > 1) {
      colnames(results)[(ncol.i+1):(ncol.i+iCols)] <-  
        paste(colnames(result.i), name, sep = "_")
    }
    if (i != m) ncol.i <- ncol.i + iCols
  }
  if (iCols == 1) {
    colnames(results) <- paste(paste(probs*100, "%", sep = ""), 
                               x$names.y, sep = "_")
  }
  return(results)
}