\name{plot}
\alias{plot.malp}
\alias{plot}
\title{Plot method for objects of class \code{"malp"}}
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
  \item{...}{Other graphical parameters}
}

\value{
  No value is returned. This method is called for producing a
  plot. It invisibly returns \code{NULL}.
}

\references{

Kim, T. & Chausse, P & Matteo Bottai & Doros, G. & Giurcanu,
M. & Luta, G. & Pena, E.A. (2026).
\dQuote{Maximum Agreement Linear Predictors.}
\emph{Electronic Journal of Statistics}, \bold{20}(2), 2892--2942.
\doi{10.1214/26-EJS2550}
  
}


\examples{

## Eye data example (see reference)

data(eye)
fitL <- malp(Stratus~Cirrus, data=subset(eye, Eye=="OS"))
plot(fitL, which="Both")

## Body fat (see reference)

data(bodyFat)
fit <- malp(PBF~ABD+WGT+FA+WRT+Age+TGH+NCK+Hip, bodyFat)
plot(fit, which="Both")


}
