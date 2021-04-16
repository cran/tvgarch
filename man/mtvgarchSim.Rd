\name{mtvgarchSim}
\alias{mtvgarchSim}

\title{Simulate from a multivariate TV-GARCH-X model}

\description{
Simulate from a multivariate multiplicative TV(s)-GARCH(p,q,r)-X model.
}

\usage{
mtvgarchSim(n, m = 2, order.g = c(1,1), order.h = c(1,1,0, 1,1,0), order.x = NULL, 
            intercept.g = c(1.2,1), size = c(0.8,1), speed = c(1.5,2), 
            location = c(0.5,0.8), intercept.h = c(0.2,0.3), arch = c(0.10,0.05), 
            garch = c(0.80,0.90), asym = NULL, xtv = NULL, xreg = NULL, par.xreg = NULL, 
            R = c(1,0.6,0.6,1), dcc = FALSE, par.dcc = NULL, 
            opt = 2, verbose = FALSE, innovations = NULL)
}

\arguments{
  \item{n}{integer}
  \item{m}{integer indicating the dimenson of the multivariate series}
  \item{order.g}{integer matrix with each row indicating the number of locations in each transition function of the TV components}
  \item{order.h}{integer matrix with each row indicating the order.h for each series; the first column controls the GARCH order, the second the ARCH order and the third the asymmetry order of the GARCH-X components}
  \item{order.x}{\code{NULL} or binary matrix indicating which xreg variables should be included as covariates in the GARCH-X components. If provided and \code{xreg} is \code{NULL}, then the selected volatility spillovers are included as covariates.}
  \item{intercept.g}{\code{NULL} or numeric vector}
  \item{size}{\code{NULL} or numeric matrix containing the size initial coefficients} 
  \item{speed}{\code{NULL} or numeric matrix containing the speed initial coefficients} 
  \item{location}{\code{NULL} or numeric matrix containing the location initial coefficients}
  \item{intercept.h}{numeric vector}
  \item{arch}{\code{NULL} or numeric matrix containing the ARCH initial coefficients}
  \item{garch}{\code{NULL} or numeric matrix containing the GARCH-type initial coefficients}
  \item{asym}{\code{NULL} or numeric matrix containing the asymmetry-type initial coefficients}
  \item{xtv}{\code{NULL} or numeric vector, time series or zoo object to include as the transition variable in the TV component. If \code{NULL}, calendar time, scaled between 0 and 1, is used as the transition variable}
  \item{xreg}{numeric vector, time series or zoo object to include as covariates in the GARCH-X component}
  \item{par.xreg}{\code{NULL} or numeric matrix containing the covariates initial coefficients}
  \item{R}{matrix of (constant) conditional correlations}
  \item{dcc}{logical. If \code{TRUE}, dynamic conditional correlations are estimated. If \code{FALSE} (default), then the conditional correlations are constant}
  \item{par.dcc}{numeric vector containing the ARCH- and GARCH-type coefficients in the dynamic conditional correlations}
  \item{opt}{integer indicating whether the speed parameter in the TV component should be scaled. If 0, no scaling; if 1, speed/sd(xtv); if 2, exp(speed)}
  \item{verbose}{\code{logical}, if \code{TRUE}, the conditional variance and innovations are also returned}
  \item{innovations}{\code{NULL} or numeric vector with the innovations. If \code{NULL}, then standard normal innovations are generated with \code{\link{rnorm}}}
}

\value{A numeric vector or matrix with the simulated values.}

\author{Susana Campos-Martins}

\seealso{
\code{\link{mtvgarch}},
\code{\link{tvgarch}},
\code{\link{garchx}},
\code{\link{zoo}}
}

\examples{
## Simulate from a bivariate CCC-TV(1)-GARCH(1,1) model (default):
mySim1 <- mtvgarchSim(n = 1500)

## Simulate from a bivariate CCC-TV(1)-GARCH(1,1)-X model (with volatility spillovers)
mySim2 <- mtvgarchSim(n = 1500, order.x = c(0,1,1,0), par.xreg = c(0,0.03,0.04,0))
}

\keyword{Time Varying Parameter Models}
\keyword{Nonlinear Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
