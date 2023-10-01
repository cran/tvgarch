\name{tvgarch-package}
\alias{tvgarch-package}
\docType{package}

\title{Time Varying GARCH Modelling}

\description{
Simulation, estimation and inference for univariate and multivariate 
TV(s)-GARCH(p,q,r)-X models, where s indicates the number and shape of the 
transition functions, p is the ARCH order, q is the GARCH order, r is the 
asymmetry order, and 'X' indicates that covariates can be included. The TV 
long-term component, as in the multiplicative TV-GARCH model of Amado and 
Ter{\"a}svirta (2013) <doi:10.1016/j.jeconom.2013.03.006>, introduces 
non-stationarity whereas the GARCH-X short-term component describes 
conditional heteroscedasticity. Maximisation by parts leads to consistent 
and asymptotically normal estimates. In the multivariate case, conditional 
variances are estimated equation by equation and dynamic conditional 
correlations are allowed.
}

\details{
  \tabular{ll}{
  Package: \tab tvgarch\cr
  Type: \tab Package\cr
  Version: \tab 2.4.1\cr
  Date: \tab 2023-10-01 \cr
  License: \tab GPL>=2\cr
  }
}

\author{
Susana Campos-Martins, \url{https://sites.google.com/site/susanacamposmartins}
\cr\cr
Maintainer: Susana Campos-Martins \cr
Contributor: Genaro Sucarrat
}

\seealso{
\code{\link{tvgarchTest}},
\code{\link{tvgarch}}, 
\code{\link{mtvgarch}}, 
\code{\link{tvgarchSim}}, 
\code{\link{mtvgarchSim}} 
}

\examples{
set.seed(123)

## Simulate from a TV(1)-GARCH(1,1) model (default):
ySim <- tvgarchSim(n = 1500)

## Test a GARCH(1,1) model against a TV(1)-GARCH(1,1) model:
yTest <- tvgarchTest(y = ySim)
yTest

## Estimate a TV(1)-GARCH(1,1) model (default):
yEst <- tvgarch(y = ySim)
yEst
}

\keyword{Time Varying Parameter Models}
\keyword{Nonlinear Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
