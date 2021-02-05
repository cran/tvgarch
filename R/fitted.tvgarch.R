fitted.tvgarch <- function (object, spec = NULL, as.zoo = TRUE, ...)
{
  if (!is.null(object$order.g)) {
    if (as.zoo == TRUE) {
      object$sigma2 <- zoo(object$sigma2, order.by = index(object$y))
      if (!is.null(object$order.g) && !is.null(spec)) {
        if (spec == "tv") object$g <- zoo(object$g, order.by = index(object$y))
        if (spec == "garch") object$h <- zoo(object$h, order.by = index(object$y))
      }
    }
    if (is.null(object$order.g)) return(object$sigma2)
    if (!is.null(object$order.g)) {
      if (is.null(spec)) return(object$sigma2)
      if (spec == "tv") return(object$g)
      if (spec == "garch") return(object$h)
    }
  }
  if (!is.null(object$order.g)) return(fitted.garchx(object = object, as.zoo = as.zoo))
}