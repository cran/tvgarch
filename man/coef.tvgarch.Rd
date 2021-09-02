\name{coef.tvgarch}

\alias{coef.tvgarch}
\alias{fitted.tvgarch}
\alias{logLik.tvgarch}
\alias{nobs.tvgarch}
\alias{predict.tvgarch}
\alias{print.tvgarch}
\alias{quantile.tvgarch}
\alias{residuals.tvgarch}
\alias{toLatex.tvgarch}
\alias{vcov.tvgarch}

\title{Extraction functions for univarate 'tvgarch' objects}

\description{Extraction functions for objects of class 'tvgarch'.}

\usage{
  \method{coef}{tvgarch}(object, spec = NULL, ...)
  \method{fitted}{tvgarch}(object, spec = NULL, as.zoo = TRUE, ...)
  \method{logLik}{tvgarch}(object, ...)
  \method{nobs}{tvgarch}(object, ...)
  \method{predict}{tvgarch}(object, n.ahead = 10, newxtv = NULL, newxreg = NULL, 
                            newindex = NULL, n.sim = 5000, verbose = FALSE, ...)
  \method{print}{tvgarch}(x, ...)
  \method{quantile}{tvgarch}(x, probs = 0.025, names = TRUE, type = 7, as.zoo = TRUE, ...)
  \method{residuals}{tvgarch}(object, as.zoo = TRUE, ...)
  \method{toLatex}{tvgarch}(object, digits = 4, ...)
  \method{vcov}{tvgarch}(object, spec = NULL, ...)
}

\arguments{
  \item{object}{an object of class 'tvgarch'.}
  \item{spec}{specifies whether the function should extract specific results. If "tv", extracts results from the TV (long-term) component; if "garch" extracts results from the GARCH-X (short-term). If \code{NULL}, extracts results from full model.}
  \item{x}{an object of class 'tvgarch'.}
  \item{as.zoo}{\code{logical}. If \code{TRUE}, then the returned result is of class \code{\link{zoo}}.}
  \item{n.ahead}{\code{integer} that determines how many steps ahead predictions should be generated.}
  \item{newxtv}{\code{NULL} or \code{vector} with the out-of-sample transition variable. If \code{NULL}, out-of-sample component g equals \code{intercept.g}.}
  \item{newxreg}{\code{vector} or \code{matrix} with the out-of-sample regressor values.}
  \item{newindex}{\code{zoo}-index for the out-of-sample predictions. If \code{NULL} (default), then \code{1:n.ahead} is used}
  \item{n.sim}{\code{integer}, the number of simulations.}
  \item{verbose}{\code{logical}. If \code{TRUE}, then the simulations - in addition to the predictions - are returned.}
  \item{probs}{\code{vector} of probabilities.}
  \item{names}{\code{logical}, whether to return names or not.}
  \item{type}{\code{integer} that determines the algorithm used to compute the quantile, see \code{\link{quantile}}.}
  \item{digits}{\code{integer}, the number of digits in the printed LaTeX code.}
  \item{\dots}{additional arguments.}
}

\value{
  \item{coef:}{numeric vector containing parameter estimates.}
  \item{fitted:}{fitted conditional variance.}
  \item{logLik:}{log-likelihood (normal density).}
  \item{nobs:}{the number of observations used in the estimation.}
  \item{predict:}{a \code{vector} with the predictions (\code{verbose=FALSE}), or a \code{matrix} with both the predictions and the simulations (\code{verbose=TRUE}).}
  \item{print:}{print of the estimation results.}
  \item{quantile:}{the fitted quantiles, i.e. the conditional standard deviation times the empirical quantile of the standardised innovations.}
  \item{residuals:}{standardised residuals.}
  \item{vcov:}{coefficient variance-covariance matrix.}
}

\references{
Cristina Amado and Timo Ter{\"a}svirta (2013) Modelling volatility by variance decomposition, Journal of Econometrics 175, 142-153. 
Cristina Amado and Timo Ter{\"a}svirta (2014) Modelling changes in the unconditional variance of long stock return series, Journal of Empirical Finance 25, 15-35. 
}

\author{Susana Campos-Martins}

\seealso{\code{\link{tvgarch}}, \code{\link{tvgarchSim}}, \code{\link{zoo}}}

\examples{
set.seed(123)

## Simulate from a TV(1)-GARCH(1,1) model (default):
ySim <- tvgarchSim(n = 1500)

## Estimate a TV(1)-GARCH(1,1) model:
yEst <- tvgarch(y = ySim)

## Print estimation results:
print(yEst)

## Extract coefficients:
coef(yEst)

## Extract and store conditional variances:
sigma2Est <- fitted(yEst)

## Extract log-likelihood:
logLik(yEst)

## Extract and store standardised residuals:
etaEst <- residuals(yEst)

## Extract variance-covariance matrix:
vcov(yEst)

## Generate predictions:
predict(yEst)
}

\keyword{Time Varying Parameter Models}
\keyword{Nonlinear Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
