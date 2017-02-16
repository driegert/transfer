#' Estimate the frequency-domain transfer function
#' 
#' Estimates a transfer function between the columns of input \code{x} and the response
#' \code{y}.
#' 
#' @param x A \code{data.frame} whose columns are the time domain input series.
#' @param y A single column \code{data.frame} containing the response series.
#' @param blockSize A \code{numeric} indicating the block sizes into which the 
#' input and response series will be partitioned.
#' @param overlap A \code{numeric} between 0 and 1 indicating how much overlap should exist 
#' between adjacent blocks.
#' @param deltat A \code{numeric} indicating the sample rate.
#' @param nw A \code{numeric} indicating the time bandwidth parameter for estimating the 
#' Slepian data tapers.
#' @param k A \code{numeric} indicating the number of tapers to use - should be approximately
#' floor(2*nw - 1) and no larger than floor(2*nw).
#' @param nFFT A \code{numeric} indicating the number of frequency bins to use (i.e. setting 
#' the zeropadding amount).
#' @param freqRange NOT CURRENTLY IMPLEMENTED.
#' @param freqOffset NOT CURRENTLY IMPLEMENTED (don't chage this... ).
#' 
#' @details Takes the times series inputs and response, divides these series into 
#' (optionally) overlapping blocks, tapers each block with Discrete 
#' Prolate Spheriodal Sequences (DPSS's or Slepian sequences), Fourier transforms each 
#' block, and then estimates the transfer function at each frequency between the Fourier 
#' transforms of the inputs and the response.
#' 
#' @return A complex matrix whose columns are the individual transfer function for each input.
#' 
#' @export
tf <- function(x, y, blockSize = dim(x)[1], overlap = 0, deltat = 1, nw = 4, k = 7, nFFT = NULL
               , freqRange = NULL, freqOffset = NULL){
  
  x2 <- sectionData(x, blockSize = blockSize, overlap = overlap)
  y2 <- sectionData(y, blockSize = blockSize, overlap = overlap)
  x3 <- taper(x2, nw = nw, k = k)
  y3 <- unlist(taper(y2, nw = nw, k = k), recursive = FALSE)
  
  time <- seq(0, (blockSize - 1)*deltat, by = deltat)
  
  if (is.null(nFFT)){
    nFFT <- 2^(floor(log2(blockSize)) + 3)
  }
  
  if (is.null(freqOffset)){
    freqOffset = 0
  }
  
  freq <- seq(0, 1/(2*deltat), by = 1/(nFFT*deltat))
  freqIdx <- 1:length(freq)
  
  H.mat <- olsTf(x = x3, y = y3, time = time, n = blockSize
             , npredictor = length(x3[[1]]), ntaper = k
             , freq = freq[freqIdx], fOffset = freqOffset)
  
  H <- data.frame(freq, H.mat)
  colnames(H) <- c("freq", colnames(x))
  
  attr(H, "class") <- c(attr(H, "class"), "transfer")
  attr(H, "blockSize") <- attr(x3, "blockSize")
  attr(H, "overlap") <- attr(x3, "overlap")
  attr(H, "numSections") <- attr(x3, "numSections")
  attr(H, "n") <- attr(x3, "n")
  attr(H, "dataTaper") <- attr(x3, "dataTaper")
  attr(H, "nw") <- attr(x3, "nw")
  attr(H, "k") <- attr(x3, "k")
  
  H
}