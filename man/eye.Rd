\name{eye}
\alias{eye}
\docType{data}
\title{Data on Optical Coherence Tomography (OCT)
}
\description{
Eye measurements used in Abedi etal. (2011). It contains two type os
measurement for 46 individuals (1 measurement per eye)
}
\usage{data("eye")}
\format{
  A data frame with 92 observations on the following 5 variables. 
  \describe{
    \item{\code{Eye}}{OD for right eye and OS for left eye.}
    \item{\code{Stratus}}{time-domain Stratus OCT.}
    \item{\code{SSS}}{Signal strength for the Stratus approach.}
    \item{\code{Cirrus}}{spectral-domain Cirrus OCT}
    \item{\code{SSC}}{Signal strength for the Cirrus approach.}    
  }
}
\references{
Abedi, G. and Patal, P. and Doros, G. and Subramanian,
M.L. (2011). Transitioning from stratus OCT to cirrus OCT: a comparison
and a proposed equation to convert central subfield macular thickness
measurements in healthy subjects. \emph{Graefe's Archive for Clinical
  and Experimental Ophthalmology}, \bold{249}, 1353--1357.
  }
\keyword{datasets}
