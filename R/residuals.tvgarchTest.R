residuals.tvgarchTest <- function (object, as.zoo = TRUE, ...)
{
  return(residuals.tvgarch(object = object$garch11, as.zoo = as.zoo))
}
