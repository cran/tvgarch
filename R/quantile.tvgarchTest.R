quantile.tvgarchTest <- function (x, probs = 0.025, names = TRUE, type = 7,
                                  as.zoo = TRUE, ...)
{
  return(quantile.tvgarch(x = x$garch11, probs = probs, names = names, 
                          type = type, as.zoo = as.zoo))
}