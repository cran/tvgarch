fitted.mtvgarch <- function (object, spec = c("sigma2", "tv", "garch", "cc"), 
                             as.zoo = TRUE, ...)
{
  spec <- match.arg(spec)
  if (spec == "sigma2") {
    if (as.zoo == TRUE) {
      object$sigma2 <- zoo(object$sigma2, order.by = index(object$y))
    }
    results <- object$sigma2
  }
  else if (spec == "tv") {
    if (is.null(object$order.g) || any(object$order.g == 0)) {
      stop ("Only GARCH models have been estimated.")
    }
    if (as.zoo == TRUE) {
      object$g <- zoo(object$g, order.by = index(object$y))
    }
    results <- object$g
  }
  else if (spec == "garch") {
    if (is.null(object$order.g) || any(object$order.g == 0)) {
      warning ("Only GARCH models have been estimated.")
    }
    if (as.zoo == TRUE) {
      object$h <- zoo(object$h, order.by = index(object$y))
    }
    results <- object$h
  }  
  else {
    if (!is.null(object$par.dcc)) {
      if (as.zoo == TRUE) {
        object$dcc <- zoo(object$dcc, order.by = index(object$y))
      }
      results = object$dcc
    }
    if (is.null(object$par.dcc)) {
      if (as.zoo == TRUE) {
        object$ccc <- zoo(object$ccc, order.by = index(object$y))
      }
      results = object$ccc
    }
  }
  return(results)
}