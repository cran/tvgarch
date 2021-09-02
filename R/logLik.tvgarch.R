logLik.tvgarch <- function(object, ...)
{
  if (!is.null(object$order.g)) {
    result <- sum(dnorm(x = object$y, mean = 0, sd = sqrt(object$h*object$g), log = TRUE))
    attr(result, "df") <- length(coef.tvgarch(object = object))
    attr(result, "nobs") <- nobs.tvgarch(object = object)
    return(result)
  }
  if (is.null(object$order.g)) return(logLik.garchx(object = object))
}
