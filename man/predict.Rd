\name{predict}
\alias{predict.malp}
\title{Predict method for object of class malp}
\description{
The function predicts the dependent variable using models estimated by
the maximum agreement method.
}
\usage{
\method{predict}{malp}(object, newdata=NULL, se.fit = FALSE,
                         interval = c("none","confidence"), level=0.95, 
                         includeLS=FALSE, LSdfCorr=FALSE,
                         vcovMet=c("Asymptotic", "Boot", "Jacknife"), ...)
}
\arguments{
  \item{object}{An object of class \code{malp}}
  \item{newdata}{An optional \code{data.frame} that contains new
    covariates used for the prediction. If not provided, the prediction
    is performed using the data used for the estimation.}
  \item{se.fit}{Should the function return the standard errors?}
  \item{interval}{If set to "confidence", a confidence interval is
    return}
  \item{level}{The level os the confidence interval.}
  \item{includeLS}{Should the method include the least squares
    prediction?}
  \item{LSdfCorr}{Should we correct for the degrees of freedom?}
  \item{vcovMet}{Method to compute the standard error of the prediction}
\item{...}{Argument for other type of objects}
}
