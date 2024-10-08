\name{vcov}
\alias{vcov.malp}
\title{Method for objects of class malp}
\description{
The function returns the covariance matrix of the coefficients obtained
using the maximum agreement method. 
}
\usage{
\method{vcov}{malp}(object, method=c("Asymptotic", "Boot", "Jackknife"),
                    B=400, LSdfCorr=FALSE, ...)
}
\arguments{
  \item{object}{An object of class \code{malp}}
  \item{method}{Method to compute the covariance matrix of the
    coefficients. By default, it is based on the asymptotic distribution
    under the normality of the data.}
  \item{LSdfCorr}{Should be apply the least squares degrees of freedom
  correction?} 
  \item{B}{The number of replications for the bootstrap method}
\item{...}{Argument for other type of objects}
}