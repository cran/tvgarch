\name{garchObj}

\alias{tv}
\alias{tvObj}
\alias{garchObj}

\title{Auxiliary functions}

\description{
Auxiliary functions used in the estimation of the univariate TV(s)-GARCH(p,q,r)-X model. Not intended for the average user.
}

\usage{
  tv(speed, location, xtv, opt, order.g)
  tvObj(par.g, fixed.par.g, xtv, opt, order.g, fixed.h, y, iter0, flag)
  garchObj(par.h, xreg, order.h, fixed.g, y, flag)
}

\arguments{
  \item{speed}{\code{NULL} or numeric vector with the values of the speed coefficients}
  \item{location}{\code{NULL} or numeric vector with the values of the location coefficients}
  \item{xtv}{\code{NULL} or numeric vector, time-series or zoo object to include as the transition variable in the TV component. If \code{NULL}, calendar time, scaled between 0 and 1, is used as the transition variable}
  \item{opt}{integer, indicates whether the speed parameter in the TV component should be scaled. If 0, no scaling; if 1, speed/sd(xtv); if 2, exp(speed).}
  \item{order.g}{integer vector of length s indicating the number of locations in each transition function of the TV component}
  \item{par.g}{numeric vector with the values of the parameters in the TV component. If \code{iter0}=TRUE, \code{par.g} takes the form c(intercept.g,size,speed,location); if \code{iter0}=FALSE, then \code{par.g}=c(size,location) and the values of the fixed parameters are provided using \code{fixed.par.g})}
  \item{fixed.par.g}{\code{NULL} or numeric vector with the values of the parameters fixed in the TV component of form c(intercept.g,speed)}
  \item{fixed.h}{numeric vector, time-series or \code{\link{zoo}} containing the values of GARCH-X component)}
  \item{y}{numeric vector, time-series or \code{\link{zoo}} object}
  \item{iter0}{logical. If \code{FALSE}, some parameters in the TV component are fixed during the iterative estimation}
  \item{flag}{integer. If 0, returns a numeric vector with the values of the objective function; if 1 returns the the value of the objective function; if 2, returns the fitted variance components}
  \item{par.h}{numeric vector with the values of the parameters in the GARCH-X component}
  \item{order.h}{integer vector of the form c(p,q,r). The first entry controls the GARCH order, the second the ARCH order and the third the asymmetry order of the GARCH-X component}
  \item{xreg}{numeric vector, time-series or zoo object to include as covariates in the GARCH-X component}
  \item{fixed.g}{numeric vector, time-series or \code{\link{zoo}} containing the values of TV component)}
}

\value{
The values of the objective function or fitted variance components
}

\author{Susana Campos-Martins}

\seealso{
\code{\link{tvgarch}}, \code{\link{fitted.tvgarch}}, \code{\link{residuals.tvgarch}} 
}

\keyword{Time-Varying Parameter Models}
\keyword{Nonlinear Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
