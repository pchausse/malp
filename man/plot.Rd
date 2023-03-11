\name{plot}
\alias{plot.malp}
\title{Plot method for objects if class malp}
\description{
The function plots the dependent variable against the fitted values from
the maximum agreement method.
}

\usage{
\method{plot}{malp}(x, y=NULL, includeLS=FALSE,
                       pch=21:22, col=2:3, bg=2:3, ...)
}
\arguments{
  \item{x}{An object of class \code{malp}}
  \item{y}{Not used at the moment.}
  \item{includeLS}{Should the method include the least squares
    prediction?}
  \item{pch}{The type of points for the MALP and LS predictions.}
  \item{col}{The colour of points for the MALP and LS predictions.}
  \item{bg}{The bachground colour of points for the MALP and LS
    predictions.}
  \item{...}{Other graphics parameters}
}
