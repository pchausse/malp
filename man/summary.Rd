\name{summary}
\alias{summary.malp}
\title{Summary for objects of class malp}
\description{
The method builds the coefficient matrices. 
}
\usage{
\method{summary}{malp}(object, vcovMet=c("Boot", "Jackknife"),
                       se=TRUE, ...)
}
\arguments{
  \item{object}{An object of class \code{malp}}
  \item{vcovMet}{Method to compute the standard error of the
    coefficients}
  \item{se}{Should the standard error of the coefficients be computed?
    By default, it is computed by simulation methods.}
\item{...}{Argument to pass to \code{vcov}.}
}
