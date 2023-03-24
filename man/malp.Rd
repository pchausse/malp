\name{malp}
\alias{malp}
\title{
Maximum Agreement Linear Predictor
}
\description{
This is the main function to estimate a linear model using the maximum
agreement predictor.
}
\usage{
malp(formula, data, obj=TRUE)
}
\arguments{
  \item{formula}{A formula for the linear model}
  \item{data}{A \code{data.frame} containing all variables included in
    the \code{formula}}
  \item{obj}{If set to \code{FALSE}, the function only returns the
    coefficients. This is faster for computing bootstrap standard
    errors} 
}


\keyword{agreement, predictor}
