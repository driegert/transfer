# testing
load("~/school_lab/contracts/health_canada2017/assets/data/edaData.RData")
x <- rnorm(1000)
y <- 2*x # need to modify this so that there is actually a relationship (needs to be done in the frequency domain?)

H <- tf(x, y, blockSize = 750, overlap = 0.25, deltat = 24*3600, nw = 6, k = 10)