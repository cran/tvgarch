logLik.tvgarch <- function(object, ...)
{
  if (!is.null(object$order.g)) {
    result <- object$logLik
    attr(result, "df") <- length(coef.tvgarch(object = object))
    attr(result, "nobs") <- nobs.tvgarch(object = object)
    return(result)
  }
  if (is.null(object$order.g)) return(logLik.garchx(object = object))
}
