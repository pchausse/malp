\name{ccc}
\alias{ccc}
\title{Concordance correlation coefficient}
\description{
This is a measure of agreement, which evaluates how close the linear
relationship between two random variables is from the 45 degree line. 
}
\usage{
ccc(x, y, type="Lin")
}
\arguments{
  \item{x}{A vector of observations.}
  \item{y}{A vector of observations with the same length as \code{x}}
  \item{type}{The type of coefficient. The only option for now is the
    coefficient intruduced by Lin (1989)}
}

\references{
Lin, L. (1989). A concordance correlation coefficient to evaluate
reproducibility. \emph{Biometrics}, \bold{45}, 255--268.
}
