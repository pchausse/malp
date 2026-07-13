\name{kFoldMalp}
\alias{kFoldMalp}
\title{
K-fold cross-validation for Maximum Agreement Linear Predictor
}
\description{
It returns some goodness of fit measures using the k-fold
cross-validation method.

}
\usage{

kFoldMalp(obj, k=5, repeated=1)

}

\arguments{
  \item{obj}{An object of class \code{"malp"}}
  \item{k}{The number of folds.}
  \item{repeated}{The number of times the k-fold method is repeated. If
  it is greater than 1, the function returns the average accuracy measures.} 
}

\value{

A numeric named vector with three goodness of fit measures: the Pearson
correlation (PCC) and the concordance correlation (CCC) coefficients between the
predictions and observed observations, and the mean squared error (MSE).

}

\keyword{k-fold cross-validation}


\references{

Kim, T. & Chausse, P & Matteo Bottai & Doros, G. & Giurcanu,
M. & Luta, G. & Pena, E.A. (2026).
\dQuote{Maximum Agreement Linear Predictors.}
\emph{Electronic Journal of Statistics}, \bold{20}(2), 2892--2942.
\doi{10.1214/26-EJS2550}
  
}

\examples{

## Stratus VS Cirrus OCT (optimal coherence tomography)

# Left eye
data(eye)
fitL <- malp(Stratus~Cirrus, data=subset(eye, Eye=="OS"))
kFoldMalp(fitL)

# Right eye
fitR <- malp(Stratus~Cirrus, data=subset(eye, Eye=="OD"))
kFoldMalp(fitR)

}
