quantile.mtvgarch <- function (x, probs = 0.025, names = TRUE, type = 7, as.zoo = TRUE, ...)
{
  results <- list()
  m <- ncol(object$y)
  nobs <- nrow(object$y)
  npar <- 0
  for(i in 1:m){
    if (!is.null(object$order.g) && object$order.g[i,1] != 0) results[[colnames(object$y)[i]]] <- quantile.tvgarch(x = object$Objs[[paste("obj", i, sep = "")]], probs = probs, names = names, type = type, as.zoo = as.zoo)
    if (is.null(object$order.g[i,1]) || object$order.g[i,1] == 0){ 
      object$Objs[[paste("obj", i, sep = "")]]$maxpqrpluss1 <- 1
      results[[colnames(object$y)[i]]] <- quantile.garchx(x = object$Objs[[paste("obj", i, sep = "")]], probs = probs, names = names, type = type, as.zoo = as.zoo)
    }
  }
  if (NCOL(results[[1]]) == 1) {
    results2 <- matrix(NA, nobs, m)
    for(i in 1:m){
      results2[,i] <- results[[i]]
    }
    results <- as.matrix(results2)
    colnames(results) <- colnames(object$y)
  }
  return(results)
}