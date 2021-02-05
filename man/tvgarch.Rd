\name{tvgarch}
\alias{tvgarch}

\title{Estimate a TV-GARCH-X model}
 
\description{ 
Quasi Maximum Likelihood (ML) estimation of a univariate multiplicative TV(s)-GARCH(p,q,r)-X model, where s indicates the number and the shape of the transition functions, r is the asymmetry order, p is the ARCH order, q is the GARCH order, and 'X' indicates that covariates can be included. Any transition variable, deterministic or stochastic, can be used to drive the transitions between the variance states. The TV long-term component introduces non-stationarity in the variance process, where the GARCH-X short-term component describes conditional heteroscedasticity. Maximization by parts leads to consistent and asymptotically normal estimates.
}

\usage{
tvgarch(y, order.g = 1, order.h = c(1,1,0), xtv = NULL, xreg = NULL, 
        initial.values = list(), opt = 2, turbo = FALSE, trace = FALSE)
}

\arguments{
\item{y}{numeric vector, time-series or \code{\link{zoo}} object}
\item{order.g}{integer vector of length s indicating the number of locations in each transition function of the TV component}
\item{order.h}{integer vector of the form c(p,q,r). The first entry controls the GARCH order, the second the ARCH order and the third the asymmetry order of the GARCH-X component}
\item{initial.values}{a list containing the initial parameter values passed on to the optimisation routines (\link{constrOptim} for the TV component and \link{nlminb} for the GARCH-X component). If list(), the default, then the values are chosen automatically. TV component: \code{intercept.g} - \code{NULL} or numeric, \code{size} - \code{NULL} or numeric vector containing the size initial coefficients, \code{speed} - \code{NULL} or numeric vector containing the speed initial coefficients, \code{location} - \code{NULL} or numeric vector containing the location initial coefficients. GARCH-X component: \code{intercept.h} - numeric, \code{arch} - \code{NULL} or numeric vector containing the ARCH initial coefficients, \code{garch} - \code{NULL} or numeric vector containing the GARCH-type initial coefficients, \code{asym} - \code{NULL} or numeric vector containing the assymetry-type initial coefficients, and \code{par.xreg} - \code{NULL} or numeric vector containing the X-type initial coefficients}
\item{xtv}{\code{NULL} or numeric vector, time-series or zoo object to include as the transition variable in the TV component. If \code{NULL}, calendar time, scaled between 0 and 1, is used as the transition variable}
\item{xreg}{numeric vector, time-series or zoo object to include as covariates in the GARCH-X component}
\item{opt}{integer indicating whether the speed parameter in the TV component should be scaled. If 0, no scaling; if 1, speed/sd(xtv); if 2, exp(speed).}
\item{turbo}{logical. If \code{FALSE} (default), then the coefficient variance-covariance is computed during estimation, and the fitted values and residuals are attached to the returned object. If \code{TRUE}, then these operations are skipped, and hence estimation is faster. Note, however, that if turbo is set to \code{TRUE}, then the coefficient-covariance, fitted values and residuals can still be extracted subsequent to estimation with vcov.tvgarch, fitted.tvgarch and residuals.tvgarch, respectively}
\item{trace}{logical. If \code{TRUE} all output is printed}
}

\value{A list of class 'tvgarch'}

\references{
Cristina Amado and Timo Ter\"asvirta (2013) Modelling volatility by variance decomposition, Journal of Econometrics 175, 142-153. Cristina Amado and Timo Ter\"asvirta (2014) Modelling changes in the unconditional variance of long stock return series, Journal of Empirical Finance 25, 15-35. 
}

\author{Susana Campos-Martins}

\seealso{
\code{\link{garchx}},
\code{\link{tvgarchSim}},
\code{\link{nlminb}},
\code{\link{constrOptim}}
}

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

\keyword{Time-Varying Parameter Models}
\keyword{Nonlinear Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
