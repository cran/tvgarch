fitted.mtvgarch <- function(object, as.zoo = TRUE, ...)
{
  if (as.zoo == TRUE) object$sigma2 <- zoo(object$sigma2, order.by = index(object$y))
  results <- list(sigma2 = object$sigma2)
  if (!is.null(object$par.dcc)) {
    if (as.zoo == TRUE) object$dcc <- zoo(object$dcc, order.by = index(object$y))
    results$dcc = object$dcc
  }
  return(results)
}