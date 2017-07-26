#' Calculates the coherence between two series
#' 
#' Estimates the frequency domain coherence using the multitaper method.
#' @param x A \code{data.frame} whose columns are the time domain input series.
#' @param y A \code{numeric} vector containing the response series.
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
#' @param forward Indicates whether the forward (TRUE) or reverse (FALSE)
#' coherence should be calculated.
#' @param average An \code{integer} representing how the average across blocks 
#' should be calculated;
#' 0 - no average, return all the individual block information; 
#' 1 - average the cross and auto-spectra across blocks, then calculate the coherency
#' 2 - estimate the coherency for each block, average the coherencey across blocks
#' 3 - estimate the MSC for each block, average the MSC across blocks.
#' @param freqRange A \code{numeric} vector containing two elements with the start 
#' and end location for the band on which to estimate the coherence.
#' @param maxFreqOffset A \code{numeric} value indicating the maximum offset coherence to 
#' calculate in the specified band.
#' @param prewhiten NOT YET IMPLEMENTED
#' @param removePeriodic NOT YET IMPLEMENTED
#' @param sigCutoff NOT YET IMPLEMENTED
#' 
#' @details MSC stands for Magnitude Squared Coherence.
#' 
#' @export
coherence2 <- function(x, y = NULL, blockSize = length(x), overlap = 0, deltat = 1
                      , nw = 4, k = 7, nFFT = NULL, forward = TRUE
                      , average = 1, msc = FALSE
                      , freqRange = NULL, maxFreqOffset = NULL
                      , prewhiten = FALSE, removePeriodic = TRUE, sigCutoff = NULL)
{
  if (!any(average == 0:3)){ 
    stop("average must have an integer value between 0 and 3.")
  }
  
  if (is.null(freqRange) && !is.null(maxFreqOffset)){
    stop("maxFreqOffset implies that freqRange should be assigned.")
  }
  
  # number of frequencies bins to use (zero-pad length)
  if (is.null(nFFT)){
    nFFT <- 2^(floor(log2(blockSize)) + 3)
  }
  
  # determine the positive values of the frequencies.
  freq <- seq(0, 1/(2*deltat), by = 1/(nFFT*deltat))
  nfreq <- length(freq)
  
  # block the data (x2, y2 are a list of data.frames)
  if (is.null(y)){
    x2 <- sectionData(data.frame(x = x[, 1], y = x[, 1]), blockSize = blockSize, overlap = overlap)
  } else {
    x2 <- sectionData(data.frame(x = x[, 1], y = y[, 1]), blockSize = blockSize, overlap = overlap)
  }
  
  numSections <- attr(x2, "numSections")
  
  if (is.null(freqRange)){
    freqIdx <- NULL
  } else {
    freqIdx <- head(which(freq >= (freqRange[1] - maxFreqOffset) ), 1):tail(which(freq <= (freqRange[2] + maxFreqOffset)), 1)
  }
  
  wtEigenCoef <- blockedEigenCoef(x2, deltat = deltat, nw = nw, k = k, nFFT = nFFT
                        , numSections = numSections, adaptiveWeighting = TRUE
                        , returnWeights = FALSE, idx = freqIdx)
  # spec <- lapply(wtEigenCoef, calculateSpec, forward = forward, idx = freqIdx)
  # spec <- blockedSpec(x2, deltat = deltat, nw = nw, k = k, nFFT = nFFT
  #                     , numSections = numSections, adaptiveWeighting = TRUE, forward = TRUE
  #                     , idx = freqIdx)
  
  
  subFreq <- freq[freqIdx]
  info = list(allFreq = freq, blockSize = blockSize, overlap = overlap
              , numSections = numSections
              , deltat = deltat, nw = nw, k = k, nFFT = nFFT
              , forward = forward, average = average, msc = msc
              , freqRange = freqRange
              , freqRangeIdx = c(head(which(freq >= freqRange[1]), 1), tail(which(freq <= freqRange[2]), 1))
              , maxFreqOffset = maxFreqOffset
              , maxFreqOffsetIdx = tail(which(freq <= maxFreqOffset), 1) - 1 #-1 due to 0 frequency
              , freqIdx = freqIdx, prewhiten = prewhiten
              , removePeriodic = removePeriodic, sigCutoff = sigCutoff)
  
  subFreqIdxRange <- c(which(freqIdx == info$freqRangeIdx[1]), which(freqIdx == info$freqRangeIdx[2]))
  
  
  if (average == 0){
    info$msc <- FALSE
    warning("Not implemented in this version.")
    return(list(freq = subFreq, coh = NULL, info = info))
  } else if (average == 1) {
    Sxy.ave <- Reduce('+', lapply(spec, "[[", "Sxy")) / numSections
    Sxx.ave <- Reduce('+', lapply(spec, "[[", "Sxx")) / numSections
    Syy.ave <- Reduce('+', lapply(spec, "[[", "Syy")) / numSections
    coh <- Sxy.ave / sqrt(Sxx.ave %*% t(Syy.ave))
  } else if (average == 2) {
    coh <- Reduce('+', lapply(spec, function(obj){ obj$Sxy / sqrt(obj$Sxx %*% t(obj$Syy)) })) / numSections
  } else if (average == 3) {
    coh <- Reduce('+', lapply(spec, function(obj){ abs(obj$Sxy)^2 / (obj$Sxx %*% t(obj$Syy)) })) / numSections
    info$msc <- TRUE
    return(list(freq = subFreq, coh = coh, info = info))
  }
  
  if (msc){
    list(freq = subFreq, coh = abs(coh)^2, info = info)
  } else {
    list(freq = subFreq, coh = coh, info = info)
  }
}

# colSums(t(A) * B)
mscOffsetByFreq <- function(offsetIdx, obj, bandIdxRange){
  lapply(obj, mscOffsetByFreqHelper, offsetIdx, bandIdxLength)
  
}

mscOffsetByFreqHelper <- function(obj, offsetIdx, bandIdxLength){
  colSums(obj$x[offsetIdx:bandIdxLength] * obj$y[offsetIdx:bandIdxLength])
}

#' Calculates the coherence between two series
#' 
#' Estimates the frequency domain coherence using the multitaper method.
#' @param x A \code{data.frame} whose columns are the time domain input series.
#' @param y A \code{numeric} vector containing the response series.
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
#' @param forward Indicates whether the forward (TRUE) or reverse (FALSE)
#' coherence should be calculated.
#' @param average An \code{integer} representing how the average across blocks 
#' should be calculated;
#' 0 - no average, return all the individual block information; 
#' 1 - average the cross and auto-spectra across blocks, then calculate the coherency
#' 2 - estimate the coherency for each block, average the coherencey across blocks
#' 3 - estimate the MSC for each block, average the MSC across blocks.
#' @param freqRange A \code{numeric} vector containing two elements with the start 
#' and end location for the band on which to estimate the coherence.
#' @param maxFreqOffset A \code{numeric} value indicating the maximum offset coherence to 
#' calculate in the specified band.
#' @param prewhiten NOT YET IMPLEMENTED
#' @param removePeriodic NOT YET IMPLEMENTED
#' @param sigCutoff NOT YET IMPLEMENTED
#' 
#' @details MSC stands for Magnitude Squared Coherence.
#' 
#' @export
offset.coh <- function(x, y = NULL, blockSize = length(x), overlap = 0, deltat = 1
                      , nw = 4, k = 7, nFFT = NULL, forward = TRUE
                      , average = 1, msc = FALSE
                      , freqRange = NULL, maxFreqOffset = NULL
                      , prewhiten = FALSE, removePeriodic = TRUE, sigCutoff = NULL)
{
  if (!any(average == 0:3)){ 
    stop("average must have an integer value between 0 and 3.")
  }
  
  if (is.null(freqRange) && !is.null(maxFreqOffset)){
    stop("maxFreqOffset implies that freqRange should be assigned.")
  }
  
  # number of frequencies bins to use (zero-pad length)
  if (is.null(nFFT)){
    nFFT <- 2^(floor(log2(blockSize)) + 3)
  }
  
  # determine the positive values of the frequencies.
  freq <- seq(0, 1/(2*deltat), by = 1/(nFFT*deltat))
  nfreq <- length(freq)
  
  # block the data (x2, y2 are a list of data.frames)
  if (is.null(y)){
    x2 <- sectionData(data.frame(x = x[, 1], y = x[, 1]), blockSize = blockSize, overlap = overlap)
  } else {
    x2 <- sectionData(data.frame(x = x[, 1], y = y[, 1]), blockSize = blockSize, overlap = overlap)
  }
  
  numSections <- attr(x2, "numSections")
  
  # wtEigenCoef <- blockedEigenCoef(x2, deltat = deltat, nw = nw, k = k
  #                                 , nFFT = nFFT, numSections = numSections
  #                                 , adaptiveWeighting = TRUE, returnWeights = FALSE)
  # 
  # if (is.null(freqRange)){
  #   spec <- lapply(wtEigenCoef, calculateSpec, forward = forward)
  #   freqIdx <- NULL
  # } else {
  #   freqIdx <- head(which(freq >= (freqRange[1] - maxFreqOffset) ), 1):tail(which(freq <= (freqRange[2] + maxFreqOffset)), 1)
  #   spec <- lapply(wtEigenCoef, calculateSpec, forward = forward, idx = freqIdx)
  # }
  
  
  if (is.null(freqRange)){
    freqIdx <- NULL
  } else {
    freqIdx <- head(which(freq >= (freqRange[1] - maxFreqOffset) ), 1):tail(which(freq <= (freqRange[2] + maxFreqOffset)), 1)
  }
  
  # spec <- lapply(wtEigenCoef, calculateSpec, forward = forward, idx = freqIdx)
  spec <- blockedSpec(x2, deltat = deltat, nw = nw, k = k, nFFT = nFFT
                      , numSections = numSections, adaptiveWeighting = TRUE, forward = TRUE
                      , idx = freqIdx)
  
  
  subFreq <- freq[freqIdx]
  info = list(allFreq = freq, blockSize = blockSize, overlap = overlap
              , numSections = numSections
              , deltat = deltat, nw = nw, k = k, nFFT = nFFT
              , forward = forward, average = average, msc = msc
              , freqRange = freqRange
              , freqRangeIdx = c(head(which(freq >= freqRange[1]), 1), tail(which(freq <= freqRange[2]), 1))
              , maxFreqOffset = maxFreqOffset
              , maxFreqOffsetIdx = tail(which(freq <= maxFreqOffset), 1) - 1 #-1 due to 0 frequency
              , freqIdx = freqIdx, prewhiten = prewhiten
              , removePeriodic = removePeriodic, sigCutoff = sigCutoff)
  
  if (average == 0){
    info$msc <- FALSE
    return(list(freq = subFreq, coh = spec, info = info))
  } else if (average == 1) {
    Sxy.ave <- Reduce('+', lapply(spec, "[[", "Sxy")) / numSections
    Sxx.ave <- Reduce('+', lapply(spec, "[[", "Sxx")) / numSections
    Syy.ave <- Reduce('+', lapply(spec, "[[", "Syy")) / numSections
    coh <- Sxy.ave / sqrt(Sxx.ave %*% t(Syy.ave))
  } else if (average == 2) {
    coh <- Reduce('+', lapply(spec, function(obj){ obj$Sxy / sqrt(obj$Sxx %*% t(obj$Syy)) })) / numSections
  } else if (average == 3) {
    coh <- Reduce('+', lapply(spec, function(obj){ abs(obj$Sxy)^2 / (obj$Sxx %*% t(obj$Syy)) })) / numSections
    info$msc <- TRUE
    return(list(freq = subFreq, coh = coh, info = info))
  }
  
  if (msc){
    list(freq = subFreq, coh = abs(coh)^2, info = info)
  } else {
    list(freq = subFreq, coh = coh, info = info)
  }
}
