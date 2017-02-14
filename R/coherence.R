## All code required for coherence

# x - dataset 1
# y - dataset 2
# xySectioned - data has already been blocked and treated.
### - x and y need to be the same length right now... otherwise this will not work.
# freqOffset - +/- freqOffset will be used to calc coherency
# slideFreq - how far to slide (in frequency) for each calculation
# slideIdx - how many frequency bins to move over to the next offset
## - set to 0 with to use W as the sliding amount
## - ignored if slideFreq is not NULL
# blockSize - how large the blocks should be
# offset for adjacent blocks (proportion = 0 < offset < 1)
# adaptive - whether to use adaptive weighting or not
# jackknife - use jackknife for variance estimates
# prewhiten - Prewhiten the blocks
# removePeriodic - removes line components (frequency domain approach)
# sigCutoff - what level of F-test to ues as the cutoff in the periodic removal
# nodes - how many cores for cluster use should this have
# clusterType - SOCK will be local, MPI would be for HPCVL

#' Calculates Frequency Domain Coherence
#' 
#' Estimates and returns the Loeve spectrum given series x and y.
#' 
#' @param x A \code{vector} of values for the first dataset.
#' @param y A \code{vector} of values for the second dataset - if null, assumes ... 
#' I'm not sure actually.
#' 
#' @return Returns a matrix containing the coherency.
#' 
#' @export
coherency.mtm <- function(x, y = NULL, xySectioned = FALSE, n = NULL, forwardCoh = TRUE
                          , freqRange = NULL, freqOffset = NULL, slideFreq = NULL, slideIdx = 1
                          , nw = NULL, k = NULL, deltat = 1, nFFT = NULL, nFFT.add = 3
                          , blockSize = NULL, overlap = NULL
                          , adaptive = TRUE, jackknife = FALSE, prewhiten = FALSE
                          , removePeriodic = TRUE, sigCutoff = NULL
                          , nodes = 1, clusterType = c("SOCK", "MPI")){
  ## Basic Setup ##
  
  if (!xySectioned){
    # not sure quite how different lengths will affect things
    nx <- length(x)
    ny <- length(y)
    
    # remove means, just in case.
    x <- x - mean(x)
    y <- y - mean(y)
  } else {
    blockSize <- length(x[[1]][[1]])
    if (is.null(n)){stop("Set n when xySectioned is TRUE.")}
    nx <- n
  }
  
  if (is.null(blockSize)){
    blockSize <- min(nx, ny)
  }
  
  if (is.null(nw)){
    nw <- 4
  }
  if (is.null(k)){
    k <- 7
  }
  if (is.null(nFFT)){
    nFFT <- 2^(floor(log2(blockSize)) + nFFT.add)
  }
  
  # how many indices to slide at each step.
  if (!is.null(slideFreq)){
    slideIdx <- round(slideFreq * deltat * nFFT)
  } else if (is.null(slideFreq) && slideIdx == 0){
    slideIdx <- floor((nw / (max(nx, ny) * deltat)) * (deltat * nFFT))
  }
  
  ## Calculate the Slepians and lambda's ##
  dpssIn <- dpss(n = blockSize, nw = nw, k = k)
  slep <- dpssIn$v
  ev <- dpssIn$eigen
  
  if (removePeriodic){
    Vk <- mvfft(rbind(slep, matrix(0, ncol = k, nrow = nFFT - blockSize)))
  }
  
  ## Don't need this right now, maybe in the future? Different blocksizes for each series? ... 
  # tmpY <- dpss(ny, nw = nw, k = k)
  # slepY <- tmpY$v
  # evY <- tmpY$eigen
  
  ## Sanity checks ##
  # Warns the user if linecomponents should be removed, but no F-test significance cutoff provided.
  if (removePeriodic && is.null(sigCutoff)){
    warning("Remove periodic components requested but no cutoff significance.  Set to (1 - 1/N).")
    sigCutoff <- 1 - 1/(min(c(length(x), length(y))))
  }
  
  # No frequency range means just a "simple" coherence at zero offset between x and y.
  if (is.null(freqRange)){
    cross <- crossSpectrum(x, y, slepX = slep, evX = ev, slepY = slep, evY = ev
                           , nFFT = nFFT, nw = nw, k = k, nFFT.add = nFFT.add
                           , adaptiveWeights = adaptive
                           , removePeriodic = removePeriodic, sigCutoff = sigCutoff
                           , returnInternals = TRUE)
    
    return(cross$Sxy / sqrt(cross$Sx * cross$Sy))
  }
  
  ## Setup Frequency Domain Indexing ##
  freq <- seq(0, 1/(2*deltat), by = 1/(deltat*nFFT))
  
  freqIdxStart <- head(which(freq >= freqRange[1]), 1)
  freqIdxEnd <- tail(which(freq <= freqRange[2]), 1)
  freqRangeIdx <- freqIdxStart:freqIdxEnd
  
  # number of indices to offset by
  numOffsetIdx <- floor(freqOffset * deltat * nFFT)
  if (slideIdx == 0){
    slideTmp <- 0
  } else {
    slideTmp <- c(rev(seq(-slideIdx, -numOffsetIdx, by = -slideIdx)), 0, seq(slideIdx, numOffsetIdx, by = slideIdx))
  }

  if (freqIdxStart - numOffsetIdx < 1 || freqIdxEnd + numOffsetIdx > (nFFT/2+1)){
    stop("Check freqRange and freqOffset - not enough frequency bins to make this combination work.")
  }
  
  ## Setup Time Domain Indexing ##
  incr <- ceiling(blockSize * (1-overlap))
  sectInd <- seq(1, nx-blockSize+1, by=incr)
  numSect <- length(sectInd)
  
  # create a list of vectors of indices - each element of the list contains the indices of the 
  ## block that will be used
  blockIdxLst <- as.list(as.data.frame(mapply(":", sectInd, sectInd+blockSize-1)))
  
  # might already be sectioned, probably not though...
  if (xySectioned){
    dat <- list()
    for (i in 1:length(x[[1]])){
      dat[[i]] <- list(x = x[[1]][[i]], y = x[[2]][[i]])
    }
  } else {
    # create a list of the data to use for parrallel
    dat <- list()
    for (i in 1:length(blockIdxLst)){
      dat[[i]] <- list(x = x[blockIdxLst[[i]]], y = y[blockIdxLst[[i]]])
    }
  }
  
  range2Start <- freqIdxStart + slideTmp
  range2End <- range2Start + length(freqRangeIdx) - 1
  # yIdxRange <- as.list(as.data.frame(mapply(":", range2Start, range2End))) # remove this - use C++ code...
  
  if (nodes <= 1){
    ## do something here - no parallelization required
    cohere <- list()
    for (i in 1:length(blockIdxLst)){
      cohere[[i]] <- coherencyHelper(dat[[i]], dpssIn = dpssIn, nw = nw, k = k, nFFT = nFFT
                                     , freqRangeIdx = freqRangeIdx, range2Start = range2Start
                                     , forwardCoh = forwardCoh
                                     , prewhiten = prewhiten, adaptive = adaptive
                                     , removePeriodic = removePeriodic, sigCutoff = sigCutoff
                                     , Vk = Vk)
    }
  } else {
    if (clusterType == "SOCK"){
      cl <- snow::makeCluster(rep("localhost", nodes), type = "SOCK")
      snow::clusterSetupRNG(cl)
      cohere <- snow::clusterApply(cl, dat, coherencyHelper
                                   , dpssIn = dpssIn, nw = nw, k = k, nFFT = nFFT
                                   , freqRangeIdx = freqRangeIdx, range2Start = range2Start
                                   , forwardCoh = forwardCoh
                                   , prewhiten = prewhiten, adaptive = adaptive,
                                   , removePeriodic = removePeriodic, sigCutoff = sigCutoff
                                   , Vk = Vk)
      snow::stopCluster(cl)
    } else if (clusterType == "MPI"){
      
    }
  }
  
  # msc <- matrix(0, nrow = dim(cohere[[1]])[1], ncol = dim(cohere[[1]])[2])
  # for (i in 1:length(cohere)){
  #   msc = msc + abs(cohere[[i]])^2
  # }
  # 
  # msc2 <- apply(msc, 1, mean)
  
  list(coherency = cohere, forwardCoh = forwardCoh
       , freqRange = freqRange, freqRangeIdx = freqRangeIdx, freqOffset = freqOffset
       , freqOffsetIdx = numOffsetIdx, slideFreq = slideFreq, slideIdx = slideIdx
       , nw = nw, k = k, deltat = deltat, nFFT = nFFT, nFFT.add = nFFT.add
       , blockSize = blockSize, overlap = overlap
       , adaptive = adaptive, jackknife = jackknife, prewhiten = prewhiten
       , removePeriodic = removePeriodic, sigCutoff = sigCutoff)
}

coherencyHelper <- function(datLst, dpssIn, nw, k, nFFT, freqRangeIdx, range2Start
                            , forwardCoh
                            , prewhiten, adaptive, removePeriodic, sigCutoff, Vk){
  n <- length(datLst$x)
  
  if (removePeriodic){
    yk.x <- removeLineCoh(spec.mtm(datLst$x, nw = nw, k = k, nFFT = nFFT, Ftest = TRUE, plot = FALSE
                                   , dpssIn = dpssIn, returnInternals = TRUE)
                   , sigCutoff = sigCutoff, Vk = Vk, k = k, nFFT = nFFT)
    yk.y <- removeLineCoh(spec.mtm(datLst$y, nw = nw, k = k, nFFT = nFFT, Ftest = TRUE, plot = FALSE
                                    , dpssIn = dpssIn, returnInternals = TRUE)
                           , sigCutoff = sigCutoff, Vk = Vk, k = k, nFFT = nFFT)
  } else {
    yk.x <- spec.mtm(datLst$x, nw = nw, k = k, nFFT = nFFT, Ftest = TRUE, plot = FALSE
                     , dpssIn = dpssIn, returnInternals = TRUE)$mtm$eigenCoefs
    yk.y <- spec.mtm(datLst$y, nw = nw, k = k, nFFT = nFFT, Ftest = TRUE, plot = FALSE
                     , dpssIn = dpssIn, returnInternals = TRUE)$mtm$eigenCoefs
  }
  
  if (prewhiten){
    # implement this... for god's sakes man!
  }
  
  if (adaptive){
    dx <- adaptiveWeightsCpp(eigenSpec = abs(yk.x)^2, ev = dpssIn$eigen, k = k, nFFT = nFFT)
    dy <- adaptiveWeightsCpp(eigenSpec = abs(yk.y)^2, ev = dpssIn$eigen, k = k, nFFT = nFFT)
  }
  
  coherencyOffsetCpp(ykx = yk.x, yky = yk.y, dx = dx, dy = dy, nTaper = k
                     , nFreqRange = length(freqRangeIdx), nOffset = length(range2Start)
                     , range2Start = range2Start, freqRangeIdx = freqRangeIdx
                     , forwardCoh = forwardCoh)
  
}

# Don't think I need this - .coherency helper will do this..
# .crossSpectrum <- function(x, y, slep, ev, nFFT, nw, k, adaptive = TRUE){
#   
# }

removeLineCoh <- function(spec, sigCutoff, Vk, k, nFFT){
  fSigIdx <- findLocalFMax(spec, cutoff = sigCutoff)
  infl <- matrix(0, nrow = nFFT/2+1, ncol = k)
  
  if (length(fSigIdx) > 0){
    for (i in 1:k){
      for (j in 1:length(fSigIdx)){
        # store the 
        ### Conjugate this because I need to use the correct side of the tapered data
        infl[, i] <- infl[, i] + shift(spec$mtm$cmv[fSigIdx[j]] * Vk[, i], -fSigIdx[j]+1)[1:(nFFT/2+1)]
      }
    }
  }
  
  eigencoef <- matrix(0, nrow = nFFT/2+1, ncol = k)
  for (i in 1:k){
    eigencoef[, i] <- spec$mtm$eigenCoefs[, i] - infl[, i]
  }
  
  eigencoef
}

# cohLst : result from coherency.mtm()
### /// fix this, should be for a specific offset.. silly David.
cohPhaseByFreq <- function(coh, offsetIdx){
  phse <- matrix(0, ncol = dim(coh$coherency[[1]])[2], nrow = length(offsetIdx))
  
  for (i in 1:length(coh$coherency)){
    phse = phse + atan2(Im(coh$coherency[[i]][offsetIdx, ]), Re(coh$coherency[[i]][offsetIdx, ]))
  }
  
  phse <- phse / length(coh$coherency)
}

mscFromCoh <- function(coh, jackknife = FALSE){
  msc <- rep(0, dim(coh$coherency[[1]])[1])
  jkVar <- rep(NA, length(msc))
    
  for (i in 1:length(coh$coherency)){
    msc <- msc + apply(abs(coh$coherency[[i]])^2, 1, mean)
  }
  
  if (jackknife){
    jkloo <- matrix(0, nrow = dim(coh$coherency[[1]])[1], ncol = dim(coh$coherency[[1]])[2])
    for (i in 1:dim(coh$coherency[[1]])[2]){
      jkloo[, i] <- apply(abs(coh$coherency[[1]][, -i])^2, 1, mean)
    }
    
    jkSlashDot <- apply(jkloo, 1, mean)
    jkVar <- (dim(jkloo)[2] - 1) * apply(matrix(mapply("-", jkloo, jkSlashDot)^2, ncol = dim(jkloo)[2])
                   , 1, mean)
  }
  
  list(msc = msc / length(coh$coherency), jkVar = jkVar)
}



freqOffsetAxis <- function(nFFT, freqOffset, deltat){
  c(rev(seq(0, -freqOffset, by = -1/(deltat * nFFT))), seq(0, freqOffset, by = 1/(deltat*nFFT))[-1])
}

offsetIdxFromFreq <- function(f, deltat, nFFT){
  f * deltat * nFFT
}
