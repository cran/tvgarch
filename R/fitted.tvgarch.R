fitted.tvgarch <- function (object, spec = c("sigma2", "tv", "garch"), 
                            as.zoo = TRUE, ...)
{
  if (!is.null(object$order.g) && object$order.g[1] != 0) {
    spec <- match.arg(spec)
    if (spec == "tv") results <- object$g 
    else if (spec == "garch") results <- object$h 
    else if (spec == "sigma2") results <- object$sigma2 
  }
  if (is.null(object$order.g) || object$order.g[1] == 0) {
    results <- object$sigma2
  }
  if (as.zoo == TRUE) {
    results <- zoo(results, order.by = index(object$y))
  }
  return(results)
}