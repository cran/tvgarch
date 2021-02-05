coef.tvgarch <- function (object, spec = NULL, ...)
{
  if (!is.null(object$order.g)) {
    if (is.null(spec)) return(c(object$par.g, object$par.h))
    if (spec == "tv") return(object$par.g)
    if (spec == "garch") return(object$par.h)
  }
  if (is.null(object$order.g)) return(coef.garchx(object = object))
}
