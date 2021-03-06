% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transfer.R
\name{predict.transfer}
\alias{predict.transfer}
\title{Predict using an estimated frequency-domain transfer function}
\usage{
\method{predict}{transfer}(object, newdata, filterMod = trim, sides = 2,
  returnInternals = FALSE, ...)
}
\arguments{
\item{object}{An object of class \code{transfer}, from a call to the \code{\link{tf}} function.}

\item{newdata}{A \code{data.frame} whose columns are the time domain input series.}

\item{filterMod}{A \code{function} to be applied to the filter coefficients before convolution.
Defaults to \code{\link{trim}}.}

\item{sides}{Argument to \code{\link{filter}}: \code{sides = 1} is for filter coefficients
that apply to past values only (a causal filter); \code{sides = 2} is for filter coefficients
centered aroung lag 0.}

\item{returnInternals}{whether to return things computed during the prediction.
Currently, this returns the full vector(s) of filter coefficients as an attribute.}

\item{...}{additional arguments passed to \code{filterMod}.}
}
\value{
A \code{data.frame} with the predicted values obtained by filtering 
the input series \code{newdata}.
}
\description{
Using output from the \code{\link{tf}} function, and an input (multivariate) time series,
this function generates predicted values.
}
\details{
The transfer function estimate is used to calculate filter coefficients to be
applied in the time domain via convolution with the input series \code{newdata}.
Prior to the convolution, the filter coefficients can be modified using the
\code{filterMod} function. If \code{filterMod} produces a causal filter, ensure that
\code{sides = 1}; in any case, be sure that the output of \code{filterMod} conforms
to the requirements of \code{\link{filter}}.

The filter coefficients are put into vectors in the orientation expected by
\code{\link{filter}}. If \code{N} is the total length of a block, the causal
coefficients from \eqn{lag = 0} to \eqn{lag = (N-1)/2} are placed on
the "right", and the non-causal coefficients from \eqn{lag = -(N-1)/(2)}
to \eqn{lag = -1} are placed on the "left". This means that for a causal filter,
you would exclude filter coefficients with an index \emph{less} than N/2,
because there are an odd number of coefficients, with \eqn{lag = 0} in the middle,
and an equal number of coefficients on each side of zero.
}
