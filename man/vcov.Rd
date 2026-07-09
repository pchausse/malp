\name{vcov}
\alias{vcov.malp}
\alias{vcov}
\title{Covariance estimation method for objects of class malp}
\description{
It returns the covariance matrix of the coefficients obtained
using the maximum agreement method. 
}
\usage{
\method{vcov}{malp}(object, method=c("Asymptotic", "Normal", "Boot", "Jackknife"),
                    B=400, LSdfCorr=FALSE, ...)
}
\arguments{
  \item{object}{An object of class \code{malp}}
  \item{method}{Method to compute the covariance matrix of the
    coefficients. See details.}
  \item{LSdfCorr}{Should we apply the least squares degrees of freedom
  correction?} 
  \item{B}{The number of replications for the bootstrap method.}
\item{...}{Argument for other type of objects.}
}

\details{

Kim et al. (2026) derive expressions for the standard errors of the
predictor directly, rather than for the covariance matrix of the
estimated coefficients. The different methods available for estimating
these standard errors are described in the Details section of
\code{\link[malp]{predict}}.

The covariance matrix of the estimated coefficients is obtained by
deriving an expression from the predictor variance formula presented in
Kim et al. (2026).
  
}


\value{

It returns the variance covariance matrix of the coefficient

}

\references{

Kim, T. & Chausse, P & Matteo Bottai & Doros, G. & Giurcanu,
M. & Luta, G. & Pena, E.A. (2026).
\dQuote{Maximum Agreement Linear Predictors.}
\emph{Electronic Journal of Statistics}, \bold{20}(2), 2892--2942.
\doi{10.1214/26-EJS2550}
  
}

\keyword{agreement predictor}

\examples{

## Stratus VS Cirrus OCT (optimal coherence tomography)

# Left eye
data(eye)
fitL <- malp(Stratus~Cirrus, data=subset(eye, Eye=="OS"))
vcov(fitL)
}


