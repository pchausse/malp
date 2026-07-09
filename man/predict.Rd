\name{predict}
\alias{predict.malp}
\alias{predict}
\title{Predict method for objects of class \code{"malp"}}
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
  \item{object}{An object of class \code{"malp"}.}
  \item{newdata}{An optional \code{data.frame} that contains new
    covariates used for the prediction. If not provided, the prediction
    is performed using the data used for the estimation.}
  \item{se.fit}{Should the function return the standard errors?}
  \item{interval}{If set to "confidence", the confidence interval is
    returned and if set to "predict", it returns the prediction intervals.}
  \item{level}{The level of the confidence interval.}
  \item{includeLS}{Should the method include the least squares
    prediction?}
  \item{LSdfCorr}{Should we correct for the degrees of freedom?}
  \item{vcovMet}{Method to compute the standard error of the prediction.}
  \item{bootInterval}{If set to \code{TRUE}, a bootstrap condifence
    interval is computed.}
  \item{bootIntType}{This is the type of bootstrap confidence
    interval. See \code{\link[boot]{boot}} for details.}
  \item{Bse.}{The number of bootstrap samples to compute the standard
    errors.}
  \item{B.}{The number of bootstrap sample to compute the confidence
    intervals.}
  \item{parallel, ncpus, cl}{These arguments are passed to
    \code{\link[boot]{boot}}.}
\item{...}{Argument for other type of objects.}
}

\details{

By default, standard errors are computed using the Delta Method, which
relies on a linear approximation of the solution. This approach is
asymptotically valid under weak regularity conditions. However, its
accuracy may be poor in small samples, in which case simulation-based
methods may provide more reliable standard error estimates.

If you are willing to assume joint normality, you can set \code{vcovMet}
to \code{"Normal"}, in which case the exact expression for the standard
error of the predictor is used. This method is valid only under the
joint normality assumption.
  
}

\value{

The following applies to all types of prediction. If \code{includeLS} is
set to \code{TRUE}, the returned value is a list with components
\code{LS} (least squares) and \code{MALP} (maximum agreement), each
containing the prediction object described below. Otherwise, the
prediction object itself is returned.

If \code{interval = "none"} and \code{se.fit = FALSE}, the prediction
object is a numeric vector of predictions.

If \code{interval = "none"} and \code{se.fit = TRUE}, the prediction
object is a list with components \code{fit}, containing the numeric
vector of predictions, and \code{se.fit}, containing the corresponding
estimated standard errors.

If \code{interval} is set to either \code{"confidence"} or
\code{"predict"}, the prediction object is a numeric matrix whose first
column contains the predictions, the second the lower bound of the
confidence or prediction interval, and the third the upper bound.

Thus, for each type of prediction, the return value follows the same
conventions as \code{\link[stats]{predict.lm}}.
  
}

\references{

Kim, T. & Chausse, P & Matteo Bottai & Doros, G. & Giurcanu,
M. & Luta, G. & Pena, E.A. (2026).
\dQuote{Maximum Agreement Linear Predictors.}
\emph{Electronic Journal of Statistics}, \bold{20}(2), 2892--2942.
\doi{10.1214/26-EJS2550}
  
}


\examples{

## Stratus VS Cirrus OCT (optimal coherence tomography)

# Left eye

data(eye)
fitL <- malp(Stratus~Cirrus, data=subset(eye, Eye=="OS"))
newd <- data.frame(Cirrus=200:205)
pr <- predict(fitL, newdata=newd, interval="confidence", includeLS=TRUE)
pr

## Body fat

data(bodyFat)
fit <- malp(PBF~ABD+WGT+FA+WRT+Age+TGH+NCK+Hip, bodyFat)
newd <- data.frame(Age=20:25, ABD=85, WGT=155, FA=27, WRT=17, TGH=60,
NCK=40, Hip=100)
pr <- predict(fit, newdata=newd, se.fit=TRUE)
pr
}