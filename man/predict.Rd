\name{predict}
\alias{predict.malp}
\title{Predict method for object of class malp}
\description{
The function predicts the dependent variable using models estimated by
the maximum agreement method.
}
\usage{
\method{predict}{malp}(object, newdata = NULL, se.fit = FALSE,
                          interval = c("none", "confidence", "prediction"), level = 0.95,
                          includeLS = FALSE, LSdfCorr = FALSE, 
                          vcovMet = c("Asymptotic", "Normal", "Boot", "Jackknife"),
                          bootInterval=FALSE,
                          bootIntType=c("all", "norm", "basic",
                                        "stud", "perc", "bca"),
                          Bse.=100, B.=300,
                          parallel = c("no", "multicore", "snow"),
                          ncpus = getOption("boot.ncpus", 1L), cl =
                          NULL, ...) 
}
\arguments{
  \item{object}{An object of class \code{malp}}
  \item{newdata}{An optional \code{data.frame} that contains new
    covariates used for the prediction. If not provided, the prediction
    is performed using the data used for the estimation.}
  \item{se.fit}{Should the function return the standard errors?}
  \item{interval}{If set to "confidence", the confidence interval is
    return and if set to "predict", it returns the predicted intervals.}
  \item{level}{The level of the confidence interval.}
  \item{includeLS}{Should the method include the least squares
    prediction?}
  \item{LSdfCorr}{Should we correct for the degrees of freedom?}
  \item{vcovMet}{Method to compute the standard error of the prediction}
  \item{bootInterval}{If set to \code{TRUE}, a bootstrap condifence
    interval is computed}
  \item{bootIntType}{This is the type of bootstrap confidence interval.}
  \item{Bse.}{The number of bootstrap samples to compute the standard
    errors.}
  \item{B.}{The number of bootstrap sample to compute the confidence
    intervals.}
  \item{parallel, ncpus, cl}{These arguments are passed to
    \code{\link[boot]{boot}}.}
\item{...}{Argument for other type of objects}
}
