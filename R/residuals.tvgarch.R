residuals.tvgarch <- function(object, as.zoo = TRUE, ...)
{
  if (!is.null(object$order.g)) {
    if (as.zoo == TRUE) object$residuals <- zoo(object$residuals, order.by = index(object$y))
    return(object$residuals)
  }
  if (is.null(object$order.g)) return(residuals.garchx(object = object, as.zoo = as.zoo))
}
