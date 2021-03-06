% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recon.R
\name{offsetRecon}
\alias{offsetRecon}
\title{Reconstructs a response based on input tranfser functions}
\usage{
offsetRecon(H, newData, hTrim = 5, sides = 2, realPart = TRUE)
}
\arguments{
\item{H}{A \code{transfer} object as returned by tf().}

\item{newData}{A \code{data.frame} containing the columns with the same name as those of \code{H}.}

\item{hTrim}{An \code{integer} indicating how points to keep on either side of the center (0 lag) point.}

\item{sides}{A value of 1 (causal) or 2 (non-causal)}

\item{realPart}{A \code{logical} indicating if the real part of the impulse response should be returned
(should be TRUE for real-valued time series).}
}
\description{
Based on multiple inputs and offset frequency coherence relationships
}
\details{
This should be incorporated into predict.transfer() at some point.  But for now, ...
}
