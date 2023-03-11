\name{summary}
\alias{summary.malp}
\title{Summary for objects of class malp}
\description{
The method builds the coefficient matrices. 
}
\usage{
\method{summary}{malp}(object, vcovMet=c("Boot", "Jackknife"), ...)
}
\arguments{
  \item{object}{An object of class \code{malp}}
  \item{vcovMet}{Method to compute the standard error of the coefficients}
\item{...}{Argument to pass to \code{vcov}.}
}
