coef.mtvgarch <- function (object, spec = c("sigma2", "tv", "garch", "cc"), ...)
{
  spec <- match.arg(spec)
  if (spec == "sigma2") {
    if (!is.null(object$order.g) && any(object$order.g != 0)) {
      results <- cbind(object$par.g, object$par.h)
    }
    else {
      results <- object$par.h
    }
  }
  else if (spec == "tv") {
    if (!is.null(object$order.g) && any(object$order.g != 0)) {
      results <- object$par.g
    }
    else {
      stop ("Only GARCH models have been estimated.")
    }
  }
  else if (spec == "garch") {
    if (!is.null(object$order.g) && any(object$order.g != 0)) {
      results <- object$par.h
    }
    else {
      results <- object$par.h
      warning ("Only GARCH models have been estimated.")
    }
  }  
  else if (spec == "cc") {
    if (is.null(object$par.dcc)) {
      results = object$ccc[1,]
    }
    if (!is.null(object$par.dcc)) {
      results$Qbar <- (1 - object$par.dcc[1] - object$par.dcc[2])*object$ccc[1,]
      results = object$par.dcc
    }
  } 
  return(results)
}