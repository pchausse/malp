\name{plot}
\alias{plot.malp}
\title{Plot method for objects if class malp}
\description{
The function plots the dependent variable against the fitted values from
the maximum agreement method. The purpose is to compare the fit with the
45 degree line representing the absolute agreement.
}
\usage{
\method{plot}{malp}(x, y=NULL, which=c("MALP", "LSLP", "Both"),
                       pch=21:22, col=2:3, bg=2:3, ...)
}
\arguments{
  \item{x}{An object of class \code{malp}}
  \item{y}{Not used at the moment.}
  \item{which}{Should we plot the fitted from MALP, LSLP or both?}
  \item{pch}{The type of points for the MALP and LS predictions.}
  \item{col}{The colour of points for the MALP and LS predictions.}
  \item{bg}{The bachground colour of points for the MALP and LS
    predictions.}
  \item{...}{Other graphics parameters}
}
