% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/old-delete.R
\name{tf.svd}
\alias{tf.svd}
\title{Estimate the frequency-domain transfer function}
\usage{
tf.svd(x, y, blockSize = dim(x)[1], overlap = 0, deltat = 1, nw = 4,
  k = 7, nFFT = NULL, freqRange = NULL, freqOffset = NULL,
  standardize = TRUE, prewhiten = TRUE, removePeriodic = TRUE)
}
\arguments{
\item{x}{A \code{data.frame} whose columns are the time domain input series.}

\item{y}{A single column \code{data.frame} containing the response series.}

\item{blockSize}{A \code{numeric} indicating the block sizes into which the 
input and response series will be partitioned.}

\item{overlap}{A \code{numeric} between 0 and 1 indicating how much overlap should exist 
between adjacent blocks.}

\item{deltat}{A \code{numeric} indicating the sample rate.}

\item{nw}{A \code{numeric} indicating the time bandwidth parameter for estimating the 
Slepian data tapers.}

\item{k}{A \code{numeric} indicating the number of tapers to use - should be approximately
floor(2*nw - 1) and no larger than floor(2*nw).}

\item{nFFT}{A \code{numeric} indicating the number of frequency bins to use (i.e. setting 
the zeropadding amount).}

\item{freqRange}{NOT CURRENTLY IMPLEMENTED.}

\item{freqOffset}{NOT CURRENTLY IMPLEMENTED (don't chage this... ).}

\item{standardize}{Should the inputs and outputs be standardized to have mean = 0 and standard deviation = 1?}
}
\value{
An object of class \code{transfer}, consisting of a complex matrix whose 
columns are the individual transfer function for each input, and several attributes
describing the transfer function estimate.
}
\description{
Estimates a transfer function between the columns of input \code{x} and the response
\code{y}. **DEPRECIATED**
}
\details{
Takes the times series inputs and response, divides these series into 
(optionally) overlapping blocks, tapers each block with Discrete 
Prolate Spheriodal Sequences (DPSS's or Slepian sequences), Fourier transforms each 
block, and then estimates the transfer function at each frequency between the Fourier 
transforms of the inputs and the response.

Different from \link{tf} in that an FFT is used and additionally, the 
standard errors of the coefficient estimates are provided as well as the 
eigenvalues from the singular value decomposition.
}
