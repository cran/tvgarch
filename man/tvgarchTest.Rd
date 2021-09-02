\name{tvgarchTest}
\alias{tvgarchTest}

\title{Test of a multiplicative time-varying GARCH model}

\description{  
Compute the non-robust and robust Lagrange-Multiplier (LM-)type test statistics for examining the null hypothesis of constant long-term variance, GARCH(1,1), against the alternative of a smoothly changing long-term component, TV-GARCH(1,1).
}

\usage{
tvgarchTest(y, xtv = NULL, alpha = 0.05, trace = TRUE)
}

\arguments{
\item{y}{numeric \code{vector}, time series or \code{\link{zoo}} object.}
\item{xtv}{\code{NULL} or numeric vector, time series or zoo object to include as the transition variable in the TV component. If \code{NULL}, calendar time, scaled between 0 and 1, is used as the transition variable.}
\item{alpha}{the significance level.}
\item{trace}{\code{logical}, if \code{TRUE} all output is printed.}
}

\value{
The number of location parameters in the single transition function.
}

\references{Cristina Amado and Timo Ter{\"a}svirta (2017) Specification and testing of multiplicative time varying GARCH models with applications, Econometric Reviews 36:4, 421-446.}

\author{Susana Campos-Martins}

\seealso{
\code{\link{tvgarch}},
\code{\link{garchx}},
\code{\link{tvgarchSim}}
}

\examples{
set.seed(123)

## Simulate from a TV(1)-GARCH(1,1) model (default):
ySim <- tvgarchSim(n = 1500)

## Test of a TV(1)-GARCH(1,1) model:
orderG1 <- tvgarchTest(y = ySim)

## Estimate a TV(1)-GARCH(1,1) model:
yEst <- tvgarch(y = ySim, order.g = orderG1)
}

\keyword{Time Varying Parameter Models}
\keyword{Nonlinear Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
