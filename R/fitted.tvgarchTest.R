fitted.tvgarchTest <- function (object, as.zoo = TRUE, ...)
{
  return(fitted.tvgarch(object = object$garch11, as.zoo = as.zoo))
}