#' Determine Offset Frequencies to include in model
#' 
#' Calculates the coherence between multiple inputs and an output and determines 
#' significant offset frequencies
#' 
#' @param dat a \code{data.frame} whose columns contain the time-domain response
#' and predictors
#' @param responseName a \code{character} string giving the name of which column is the 
#' response.
#' @param blockSize The length of a single block to use (if blocking)
#' @param overlap A \code{numeric} value in the range [0, 1) indicating the proporation of 
#' overlap between neighbouring blocks.
#' @param deltat The sampling rate in seconds.
#' @param nw time-bandwidth parameter for multitaper
#' @param k number of tapers to use (k < 2*nw)
#' @param nFFT the number of frequency bins to use (nFFT > 2*ndata)
#' @param freqRange A vector with 2 elements containing the start and end frequencies (in Hz) 
#' over which to calculate the coherence.
#' @param maxFreqOffset Every pair of frequencies between f1 (series 1) 
#' and f1 +/- maxFreqOffset (series 2) will be calculated (f1 + maxFreqOffset < nyquist)
#' @param calcType An \code{integer} value indicating how averaging over blocks should 
#' be performed:
#' 1 - calculate the MSC on each block, then average;
#' 2 - calculate the cross and auto spectra on each block, average each quantity 
#' across blocks, then calculate the coherency;
#' 3 - calculate the coherency on each block, then average
#' 4 - minimum MSC across blocks
#' @param forward An \code{integer} indicating whether the forward (1) or reverse (0) coherence
#' should be calculated.
#' 
#' @export
offsetFreq <- function(dat, responseName, blockSize = dim(dat)[1], overlap = 0
                       , deltat = 1, nw = 4, k = 7, nFFT = NULL, mscCutoff = 0.2
                       , freqRange = NULL, maxFreqOffset = 0
                       , useZeroOffset = TRUE, nOffsetFreq = -1
                       , calcType = 1, forward = 1){
  if (is.null(freqRange) || maxFreqOffset == 0){
    stop("freqRange or maxFreqOffset are not set - no point in running this.")
  }
  
  pNames <- names(dat)
  if (is.na(match(responseName, pNames))){
    stop("'responseName' does not correspond to a column name of 'dat'.")
  }
  
  pNames <- pNames[-match(responseName, pNames)]
  
  if (is.null(nFFT)){
    nFFT <- 2^(floor(log2(blockSize)) + 3)
  }
  
  for (i in 1:length(pNames)){
    offCoh <- transfer2::coherence(dat[, responseName], dat[, pNames[i] ]
                                   , blockSize = blockSize
                                   , overlap = overlap, dt = deltat
                                   , nw = nw, k = k, nFFT = nFFT
                                   , freqRange = freqRange
                                   , maxFreqOffset = maxFreqOffset
                                   , calcType = calcType
                                   , name1 = responseName, name2 = pNames[i])
    if (i == 1){
      offsets <- maxOffCoh(offCoh, cutoff = mscCutoff, nw = nw, N = blockSize, nFFT = nFFT
                           , name = offCoh$info$d2, useZeroOffset = useZeroOffset
                           , nOffsetFreq = nOffsetFreq)
    } else {
      offsets <- maxOffCoh(offCoh, cutoff = mscCutoff, nw = nw, N = blockSize, nFFT = nFFT
                           , name = offCoh$info$d2, useZeroOffset = useZeroOffset
                           , nOffsetFreq = nOffsetFreq, appendTo = offsets)
    }
  }
  
  info <- list(response = responseName, predictors = pNames, ndata = dim(dat)[1]
               , blockSize = blockSize, overlap = overlap, deltat = deltat
               , nw = nw, k = k, nFFT = nFFT, freqRange = freqRange
               , freqRangeIdx = offCoh$info$freqRangeIdx
               , maxFreqOffset = maxFreqOffset
               , maxFreqOffsetIdx = offCoh$info$maxFreqOffsetIdx
               , calcType = calcType, calcTypeDesc = offCoh$info$calcTypeDesc
               , mscCutoff = mscCutoff
               , forward = forward)
  list(offsets = offsets, info = info)
}

# takes an object as returned from transfer2::coherence()
### NOTE: divide by 4 here in the main lobe width - should that be there?
maxOffCoh <- function(coh, cutoff, nw, N, nFFT, name = "d2"
                      , useZeroOffset = TRUE, nOffsetFreq = -1, appendTo = NULL){
  if (useZeroOffset){
    lobeWidthIdx <- ceiling(nFFT / N) # Rayleigh number of bins
    #lobeWidthIdx <- ceiling((nw / N) * (nFFT) / 4) # how many bins make up W - don't use anything in here close to 0 offset.
    
  } else {
    lobeWidthIdx <- 0
  }
  
  zeroOff <- ceiling(length(coh$offset)/2)
  coh.df <- as.data.frame(coh$coh)
  # standard use - includes coherences above cutoff "away" from zero-offset
  if (useZeroOffset & nOffsetFreq == -1){
    idxOff <- lapply(coh.df
                     , function(x, c){
                       n <- length(x)
                       tmp <- (which(x[2:(n-1)] > x[1:(n-2)] 
                                     & x[2:(n-1)] > x[3:n] 
                                     & x[2:(n-1)] > c) + 1) - zeroOff
                       tmp[which(!(tmp >= -lobeWidthIdx & tmp <= lobeWidthIdx))]
                     }
                     , c = cutoff)
  # only includes the nOffsetFreq highest coherence offsets.
  } else if (useZeroOffset & nOffsetFreq != -1){
    idxOff <- lapply(coh.df
                     , function(x, c){
                       n <- length(x)
                       tmp <- (which(x[2:(n-1)] > x[1:(n-2)] 
                                     & x[2:(n-1)] > x[3:n] 
                                     & x[2:(n-1)] > c) + 1) - zeroOff
                       tmp2 <- tmp[which(!(tmp >= -lobeWidthIdx & tmp <= lobeWidthIdx))]
                       tail( tmp2[order(x[tmp2 + zeroOff]) ], nOffsetFreq)
                     }
                     , c = cutoff)
  # Uses at least 1 value (max coherence), regardless of cutoff, up to a maximum of nOffsetFreq
  } else if (!useZeroOffset){
    idxOff <- lapply(coh.df
                     , function(x, c){
                       n <- length(x)
                       tmp <- (which(x[2:(n-1)] > x[1:(n-2)] 
                                     & x[2:(n-1)] > x[3:n] 
                                     & x[2:(n-1)] > c) + 1) - zeroOff
                       if (length(tmp) == 0){
                         return(which(x == max(x)) - zeroOff)
                       } else if (nOffsetFreq != -1) {
                         return(tail( tmp[order(x[tmp + zeroOff]) ], nOffsetFreq))
                       } else {
                         return(tmp)
                       }
                     }
                     , c = cutoff)
  }
  
  # names(idxOff) <- rep("offIdx", length(idxOff))
  
  freqRangeIdx <- coh$info$freqRangeIdx[1]:coh$info$freqRangeIdx[2]
  if (is.null(appendTo)){
    # create the list and sublists if need be ... 
    appendTo <- list()
    for (i in 1:length(idxOff)){
      if (length(idxOff[[i]]) == 0) {
        appendTo[[ paste0("fBin", freqRangeIdx[i]) ]] <- list()
        appendTo[[ paste0("fBin", freqRangeIdx[i]) ]][[name]] <- list(offIdx = idxOff[[i]]
                                                                      , offFreq = numeric(0)
                                                                      , offCoh = numeric(0))
      } else {
        appendTo[[ paste0("fBin", freqRangeIdx[i]) ]] <- list()
        appendTo[[ paste0("fBin", freqRangeIdx[i]) ]][[name]] <- list(offIdx = idxOff[[i]]
                                                                      , offFreq = coh$offset[ zeroOff + idxOff[[i]] ]
                                                                      , offCoh = coh$coh[zeroOff+idxOff[[i]], i ])
      }
    }
  } else {
    for (i in 1:length(idxOff)){
      if (length(idxOff[[i]]) == 0) {
        appendTo[[ paste0("fBin", freqRangeIdx[i]) ]][[name]] <- list(offIdx = idxOff[[i]]
                                                                      , offFreq = numeric(0)
                                                                      , offCoh = numeric(0))
      } else {
        appendTo[[ paste0("fBin", freqRangeIdx[i]) ]][[name]] <- list(offIdx = idxOff[[i]]
                                                                      , offFreq = coh$offset[ zeroOff + idxOff[[i]] ]
                                                                      , offCoh = coh$coh[zeroOff+idxOff[[i]], i ])
      }
    }
  }
  appendTo
}


# creates a large, but sparse, matrix containing the transfer functions with offsets.
# obj - is H.tmp from tf()
offsetTfMatrix <- function(obj){
  names <- unique(unlist(lapply(obj, function(x) colnames(x$coef))))
  
  H <- matrix(0, nrow = length(obj), ncol = length(names))
  colnames(H) <- names
  for (i in 1:length(obj)){
    H[i, colnames(obj[[i]]$coef)] <- obj[[i]]$coef[1, ]
  }
  
  H
}