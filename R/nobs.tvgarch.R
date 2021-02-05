nobs.tvgarch <- function (object, ...)
{
  if (!is.null(object$order.g)) return(length(object$sigma2))
  if (is.null(object$order.g)) return(nobs.garchx(object = object))
}
