\name{transfer-package}
\alias{transfer-package}
\alias{transfer}
\docType{package}
\title{
  Estimating frequency-domain transfer functions
}
\description{
  Estimating frequency-domain transfer functions.
}
\details{
  Frequency-domain transfer functions are linear regressions of
  a Fourier-transformed \emph{response} time series onto a set of
  Fourier-transformed \emph{input} time series. The regression is
  estimated at a set of frequencies. To obtain additional samples,
  the input and response series are tapered using multiple
  orthogonal tapers (discrete prolate spheroidal sequences), and/or
  estimated across time-blocks. The function \link{tf} is the main
  function.
}
\author{
David Riegert

Maintainer: David Riegert <driegert@gmail.com>
}
\references{
Chave, A. D.,  and Jones, A. G. 2012 The magnetotelluric method: Theory and practice.
Cambridge University Press.

Thomson, D.J. 1982 Spectrum estimation and harmonic analysis. Proceedings of the IEEE,
70:1055–1096.
}
\keyword{ package }
