\name{eye}
\alias{eye}
\docType{data}
\title{Data on Optical Coherence Tomography (OCT)
}
\description{
Eye measurements used in Abedi et al. (2011). The original dataset
contains two type of measurement for 46 individuals (1 measurement per
eye). The data is then filtered to include only observations for which
the signal strength is greater or equal to 6 for both measurement methods.
}
\usage{data("eye")}
\format{
  A data frame with 56 observations on the following 3 variables. 
  \describe{
    \item{\code{Eye}}{OD for right eye and OS for left eye.}
    \item{\code{Stratus}}{time-domain Stratus OCT.}
    \item{\code{Cirrus}}{spectral-domain Cirrus OCT}
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
