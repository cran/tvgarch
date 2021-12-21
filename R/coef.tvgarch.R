coef.tvgarch <- function (object, spec = c("sigma2", "tv", "garch"), ...)
{
  if (!is.null(object$order.g) && object$order.g[1] != 0) {
    spec <- match.arg(spec)
    if (spec == "sigma2") results <- c(object$par.g, object$par.h)
    else if (spec == "tv") results <- object$par.g
    else if (spec == "garch") results <- object$par.h
  } 
  if (is.null(object$order.g) || object$order.g[1] == 0) results <- object$par.h
  return(results)
}
