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
    coefficient introduced by Lin (1989)}
}

\references{
Lin, L. (1989). A concordance correlation coefficient to evaluate
reproducibility. \emph{Biometrics}, \bold{45}, 255--268.
}

\examples{
## High concordance: almost y=x
set.seed(112233)
x <- rnorm(100)
u <- rnorm(100, sd=0.1)
y1 <- x + u
ccc(x,y1)

## Lower concordance: almost y=1+x
y2 <- 1+x+u
ccc(x,y2)

## But same Pearson correlation
cor(x,y2)

}