\name{vcov}
\alias{vcov.malp}
\title{Method for objects of class malp}
\description{
The function returns the covariance matrix of the coefficients obtained
using the maximum agreement method. 
}
\usage{
\method{vcov}{malp}(object, method=c("Boot", "Jackknife"), B=400, ...)
}
\arguments{
  \item{object}{An object of class \code{malp}}
  \item{method}{Method to compute the covariance matrix of the
    coefficients}
  \item{B}{The number of replications for the bootstrap method}
\item{...}{Argument for other type of objects}
}