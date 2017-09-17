#' Estimate the frequency-domain transfer function
#' 
#' Estimates a transfer function between the columns of input \code{x} and the response
#' \code{y}. **DEPRECIATED**
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
#' @param standardize Should the inputs and outputs be standardized to have mean = 0 and standard deviation = 1? 
#' 
#' @details Takes the times series inputs and response, divides these series into 
#' (optionally) overlapping blocks, tapers each block with Discrete 
#' Prolate Spheriodal Sequences (DPSS's or Slepian sequences), Fourier transforms each 
#' block, and then estimates the transfer function at each frequency between the Fourier 
#' transforms of the inputs and the response.
#' 
#' Different from \link{tf} in that an FFT is used and additionally, the 
#' standard errors of the coefficient estimates are provided as well as the 
#' eigenvalues from the singular value decomposition.
#' 
#' @return An object of class \code{transfer}, consisting of a complex matrix whose 
#' columns are the individual transfer function for each input, and several attributes
#' describing the transfer function estimate.
#' 
#' @export
tf.svd <- function(x, y, blockSize = dim(x)[1], overlap = 0, deltat = 1, nw = 4
                   , k = 7, nFFT = NULL, freqRange = NULL, freqOffset = NULL
                   , standardize = TRUE, prewhiten = TRUE, removePeriodic = TRUE){
  warning("tf.svd() is a depreciated function - use tf() with method == 'svd' instead.")
  # standardize the series by removing the mean and dividing by the standard deviation
  if( standardize ){
    stdPars <- vector( mode = "list" )
    stdPars$x <- data.frame( xmean = sapply( x, mean ), xsd = sapply( x, sd ) )
    stdPars$y <- data.frame( ymean = sapply( y, mean ), ysd = sapply( y, sd ) )
    std <- function( a ) (a - mean(a))/sd(a)
    x <- data.frame( lapply( x, std ) )
    y <- data.frame( lapply( y, std ) )
  }
  
  # I *think* we want to prewhiten based on the most data available?
  # (probably.. due to welch estimate, etc... )
  
  # block the data (x2, y2 are a list of data.frames)
  x2 <- sectionData(x, blockSize = blockSize, overlap = overlap)
  y2 <- sectionData(y, blockSize = blockSize, overlap = overlap)
  
  # determine an appropriate time vector (0 to T-dt)
  time <- seq(0, (blockSize - 1)*deltat, by = deltat)
  
  # number of frequencies bins to use (zero-pad length)
  if (is.null(nFFT)){
    nFFT <- 2^(floor(log2(blockSize)) + 3)
  }
  
  # determine the positive values of the frequencies.
  freq <- seq(0, 1/(2*deltat), by = 1/(nFFT*deltat))
  
  nfreq <- length(freq)
  
  # freqIdx at which frequencies of the response the transfer function should be 
  # estimated (you could not use all frequencies if offsets were also used)
  # NOT USED in this iteration of tf() ... 
  # freqIdx <- 1:length(freq)
  
  # taper the blocked data
  # x3 <- taper(x2, nw = nw, k = k)
  # y3 <- unlist(taper(y2, nw = nw, k = k), recursive = FALSE)
  
  # I don't think we actually need to store ALL of this (at least not the spec.mtm objects..)
  x.spec <- list()
  x.wtEigenCoef <- list()
  y.spec <- list()
  y.wtEigenCoef <- list()
  
  numSections <- attr(x2, "numSections")
  
  # used in the next for-loop ... 
  weightedEigen <- function(obj){ obj$mtm$eigenCoefs * obj$mtm$eigenCoefWt }
  
  for (i in 1:numSections){
    x.spec[[i]] <- lapply(x2[[i]], spec.mtm, deltat = deltat, dtUnits = "second", nw = nw
                          , k = k, nFFT = nFFT, plot = FALSE, returnInternals = TRUE)
    x.wtEigenCoef[[i]] <- lapply(x.spec[[i]], weightedEigen)
    y.spec[[i]] <- lapply(y2[[i]], spec.mtm, deltat = deltat, dtUnits = "second", nw = nw
                          , k = k, nFFT = nFFT, plot = FALSE, returnInternals = TRUE)
    y.wtEigenCoef[[i]] <- lapply(y.spec[[i]], weightedEigen)
  }
  
  
  # indexing helper function that grabs and stacks all the eigencoefficients at a frequency
  eigenByFreq <- function(obj, rowNum, numEl){
    matrix(unlist(lapply(obj, function(x, idx) x[idx, ], rowNum)), ncol = numEl)
  }
  
  
  x.design <- list()
  # y.design <- list()
  
  # stack the eigencoefficients by frequency from each block
  # form into a single list
  for (i in 1:nfreq){
    x.design[[i]] <- list(x = do.call(rbind, lapply(x.wtEigenCoef, eigenByFreq, rowNum = i, numEl = 3))
                          , y = do.call(rbind, lapply(y.wtEigenCoef, eigenByFreq, rowNum = i, numEl = 1)))
    # y.design[[i]] <- do.call(rbind, lapply(y.wtEigenCoef, eigenByFreq, rowNum = i, numEl = 1))
  }
  
  # This would be used if incorporating offset frequencies into the transfer function
  # regression.
  # NOT IMPLEMENTED YET (and not needed for the calling function)
  # if (is.null(freqOffset)){
  #   freqOffset = 0
  # }
  
  # perform the least squares complex-regresison.
  # H.mat <- olsTf(x = x3, y = y3, time = time, n = blockSize
  #                , npredictor = length(x3[[1]]), ntaper = k
  #                , freq = freq[freqIdx], fOffset = freqOffset)
  
  H.tmp <- lapply(x.design, function(obj){ svdRegression(obj$x, obj$y) })
  
  H.mat <- do.call(rbind, lapply(H.tmp, "[[", "coef"))
  # stdErr <- do.call(rbind, lapply(H.tmp, "[[", "stdErr"))
  # svdEigenval <- do.call(rbind, lapply(H.tmp, "[[", "ev"))
  
  H <- data.frame(freq, H.mat)
  colnames(H) <- c("freq", colnames(x))
  
  
  attr(H, "class") <- c(attr(H, "class"), "transfer")
  attr(H, "stdErr") <- do.call(rbind, lapply(H.tmp, "[[", "stdErr"))
  attr(H, "svdEv") <- do.call(rbind, lapply(H.tmp, "[[", "ev"))
  attr(H, "blockSize") <- attr(x2, "blockSize")
  attr(H, "overlap") <- attr(x2, "overlap")
  attr(H, "numSections") <- attr(x2, "numSections")
  attr(H, "n") <- attr(x2, "n")
  attr(H, "dataTaper") <- "Slepian"
  attr(H, "nw") <- nw
  attr(H, "k") <- k
  attr(H, "standardize") <- standardize
  if( standardize ) attr(H, "stdPars") <- stdPars
  
  # I made the standard errors and eigenvalues attributes on the transfer function
  # so that predict method wouldn't need any altering.  I don't know if this is the 
  # best way to do it though.
  H
}