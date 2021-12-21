predict.tvgarchTest <- function (object, n.ahead = 10, newxreg = NULL, 
                                 newindex = NULL, n.sim = 5000, as.zoo = TRUE, 
                                 verbose = FALSE, ...)
{
  return(predict.tvgarch(object = object$garch11, n.ahead = n.ahead, 
                         newxreg = newxreg, newindex = newindex, n.sim = n.sim, 
                         as.zoo = as.zoo, verbose = verbose))
}