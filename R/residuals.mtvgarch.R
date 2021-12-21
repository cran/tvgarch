residuals.mtvgarch <- function (object, as.zoo = TRUE, ...)
{
  if (as.zoo == TRUE) {
    object$residuals <- zoo(object$residuals, order.by = index(object$y))
  }
  return(object$residuals)
}
