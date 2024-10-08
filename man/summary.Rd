\name{summary}
\alias{summary.malp}
\title{Summary for objects of class malp}
\description{
The method builds the coefficient matrices. 
}
\usage{
\method{summary}{malp}(object, vcovMet=c("Asymptotic","Boot", "Jackknife"),
                       se=TRUE, LSdfCorr=FALSE, ...)
}
\arguments{
  \item{object}{An object of class \code{malp}}
  \item{vcovMet}{Method to compute the covariance matrix of the
    coefficients. By default, it is based on the asymptotic distribution
    under the normality of the data.}  
  \item{se}{Should the standard error of the coefficients be computed?}
  \item{LSdfCorr}{Should be apply the least squares degrees of freedom
    correction?} 
\item{...}{Argument to pass to \code{vcov}.}
}
