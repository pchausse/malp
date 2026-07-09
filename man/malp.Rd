\name{malp}
\alias{malp}
\title{
Maximum Agreement Linear Predictor
}
\description{
This is the main function to estimate the maximum agreement linear
predictor. Instead of estimating the predictor directly, as in Kim et
al. (2026), the coefficients are estimated like linear models fitted
by \code{\link[stats]{lm}} and saved in the \code{"malp"} object created
by this function.  

}
\usage{
malp(formula, data, obj=TRUE)
}

\arguments{
  \item{formula}{A formula for the linear model}
  \item{data}{A \code{data.frame} containing all variables included in
    the \code{formula}}
  \item{obj}{If set to \code{FALSE}, the function only returns the
    coefficients. It saves a little bit of time for large models if all
  we need is the coefficients (e.g. in simulation studies). } 
}

\value{

\code{malp} returns an object of class \code{"malp"}, with the following components:
    \item{coefficients}{Coefficient estimates.}
    \item{gamma}{Correlation between the dependent variable and the
    least squares fitted values.}
    \item{lm}{Least squares object fitted by \code{\link[stats]{lm}}.}
    \item{varY}{Estimated variance of the dependent variable.}
    \item{muY}{Sample mean of the dependent variable.}
    \item{call}{The original function call.}
    \item{data}{The data set used for the estimation. The missing values
    are therefore omitted.}
    \item{na.action}{A vector indicating the rows that were omitted from
    the dataset due to missing values. It is \code{NULL}
    if the dataset does not contain missing values.}

The \code{print} method for this object is like the method for objects
of class \code{"lm"}. It prints the call and coefficient estimates.
  
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
fitL
# Right eye
fitR <- malp(Stratus~Cirrus, data=subset(eye, Eye=="OD"))
fitR

## Body fat

data(bodyFat)
fit <- malp(PBF~ABD+WGT+FA+WRT+Age+TGH+NCK+Hip, bodyFat)
fit

}