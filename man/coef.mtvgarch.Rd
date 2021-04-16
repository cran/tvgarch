\name{coef.mtvgarch}

\alias{coef.mtvgarch}
\alias{fitted.mtvgarch}
\alias{logLik.mtvgarch}
\alias{nobs.mtvgarch}
\alias{predict.mtvgarch}
\alias{print.mtvgarch}
\alias{quantile.mtvgarch}
\alias{residuals.mtvgarch}
\alias{toLatex.mtvgarch}
\alias{vcov.mtvgarch}

\title{Extraction functions for multivariate 'mtvgarch' objects}

\description{Extraction functions for objects of class 'mtvgarch'}

\usage{
  \method{coef}{mtvgarch}(object, ...)
  \method{fitted}{mtvgarch}(object, as.zoo = TRUE, ...)
  \method{logLik}{mtvgarch}(object, ...)
  \method{nobs}{mtvgarch}(object, ...)
  \method{predict}{mtvgarch}(object, n.ahead = 10, newxtv = NULL, newxreg = NULL, 
  newindex = NULL, ...)
  \method{print}{mtvgarch}(x, ...)
  \method{quantile}{mtvgarch}(x, probs = 0.025, names = TRUE, type = 7, as.zoo = TRUE, ...)
  \method{residuals}{mtvgarch}(object, as.zoo = TRUE, ...)
  \method{toLatex}{mtvgarch}(object, digits = 4, ...)
  \method{vcov}{mtvgarch}(object, ...)
}

\arguments{
  \item{object}{an object of class 'tvgarch'}
  \item{x}{an object of class 'tvgarch'}
  \item{as.zoo}{logical. If \code{TRUE}, then the returned result is of class \code{\link{zoo}}}
  \item{n.ahead}{\code{integer} that determines how many steps ahead predictions should be generated}
  \item{newxtv}{\code{NULL} or \code{vector} with the out-of-sample transition variable. If \code{NULL}, out-of-sample component g equals \code{intercept.g}}
  \item{newxreg}{\code{vector} or \code{matrix} with the out-of-sample regressor values}
  \item{newindex}{\code{zoo}-index for the out-of-sample predictions. If \code{NULL} (default), then \code{1:n.ahead} is used}
  \item{probs}{\code{vector} of probabilities}
  \item{names}{\code{logical}, whether to return names or not}
  \item{type}{\code{integer} that determines the algorithm used to compute the quantile, see \code{\link{quantile}} }
  \item{digits}{\code{integer}, the number of digits in the printed LaTeX code}
  \item{\dots}{additional arguments}
}

\value{
  \item{coef:}{list of numeric vectors containing parameter estimates}
  \item{fitted:}{fitted conditional variances and correlations}
  \item{logLik:}{list of log-likelihood values (normal density)}
  \item{nobs:}{number of observations used in the estimation}
  \item{predict:}{matrix with the variance predictions}
  \item{print:}{print of the estimation results}
  \item{quantile:}{the fitted quantiles, i.e. the conditional standard deviation times the empirical quantile of the standardised innovations}
  \item{residuals:}{standardised residuals}
  \item{vcov:}{list of coefficient variance-covariance matrices}
}

\references{
Cristina Amado and Timo Ter\"{a}svirta (2013) Modelling volatility by variance decomposition, Journal of Econometrics 175, 142-153. 
Christian Francq and Jean-Michel Zako\"{i}an (2016) Estimating multivariate volatility models equation by equation, J. R. Stat. Soc. Ser. B Stat. Methodol 78, 613-635.
}

\author{Susana Campos-Martins}

\seealso{
\code{\link{mtvgarch}}, 
\code{\link{mtvgarchSim}}, 
\code{\link{tvgarch}}, 
\code{\link{garchx}}, 
\code{\link{zoo}}
}

\examples{
set.seed(123)

## Simulate from a bivariate CCC-TV(1)-GARCH(1,1) model (default):
mySim <- mtvgarchSim(n = 1500)

## Estimate a CCC-TV(1)-GARCH(1,1) model:
myEst <- mtvgarch(y = mySim)

## Print estimation results:
print(myEst)

## Extract coefficients:
coef(myEst)

## Extract and store conditional variances:
sigma2Est <- fitted(myEst)

## Extract log-likelihood:
logLik(myEst)

## Extract and store standardised residuals:
etaEst <- residuals(myEst)

## Extract variance-covariance matrix:
vcov(myEst)

## Generate predictions:
predict(myEst)
}

\keyword{Time Varying Parameter Models}
\keyword{Nonlinear Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
