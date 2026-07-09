\name{summary}
\alias{summary}
\alias{summary.malp}
\title{Summary for objects of class \code{"malp"}}
\description{
The method computes a coefficient matrix similar to the summary method
for \code{"lm"} objects. 
}
\usage{
\method{summary}{malp}(object,
                       vcovMet=c("Asymptotic", "Normal", "Boot", "Jackknife"),
                       se=TRUE, LSdfCorr=FALSE, ...)
}
\arguments{
  \item{object}{An object of class \code{"malp"}}
  \item{vcovMet}{Method to compute the covariance matrix of the
    coefficients. See the details in \code{\link[malp]{vcov}}.}
  \item{se}{Should the standard error of the coefficients be computed?
    This argument was included when \code{"Boot"} was the default
    method. It is no longer needed, but we keep it for now.}
  \item{LSdfCorr}{Should we apply the least squares degrees of freedom
    correction?} 
\item{...}{Argument to pass to \code{\link[malp]{vcov}}.}
}

\value{

\code{symmary.malp} returns an object of class \code{"summary.malp"}, with the following components:
    \item{coefficients}{The coefficient matrix with estimates, standard
    errors, t-ratios and p-values.}
    \item{fitMALP}{Measures of goodness of fit for maximum agreement predictors.}
    \item{fitLSLP}{Measures of goodness of fit for least squares predictors.}
    \item{call}{The original function call.}

A \code{print} method is also available for this class of objects.
  
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
summary(fitL)
# Right eye
fitR <- malp(Stratus~Cirrus, data=subset(eye, Eye=="OD"))
summary(fitR)

## Body fat

data(bodyFat)
fit <- malp(PBF~ABD+WGT+FA+WRT+Age+TGH+NCK+Hip, bodyFat)
summary(fit)

}
