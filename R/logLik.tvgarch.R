logLik.tvgarch <- function (object, ...)
{
    result <- object$logLik
    class(result) <- "logLik"
    attr(result, "df") <- length(coef.tvgarch(object = object))
    attr(result, "nobs") <- length(object$sigma2)
    return(result)
}
