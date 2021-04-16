coef.mtvgarch <- function(object, ...)
{
  results <- list()
  m <- ncol(object$y)
  for(i in 1:m){
    object.i <-object$Objs[[paste("obj", i, sep = "")]]
    if (!is.null(object$order.g) && object$order.g[i,1] != 0) results[[paste("tvgarch", i, sep = "")]] <- coef.tvgarch(object = object.i)
    if (is.null(object$order.g[i,1]) || object$order.g[i,1] == 0){ 
      results[[paste("garch", i, sep = "")]] <- coef.garchx(object = object.i)
    }
  }  
  if (is.null(object$par.dcc)) results$ccc = object$ccc[1,]
  if (!is.null(object$par.dcc)) {
    results$Qbar <- (1 - object$par.dcc[1] - object$par.dcc[2]) * object$ccc[1,]
    results$dcc = object$par.dcc
  }
  return(results)
}