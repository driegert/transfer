#include <Rcpp.h>
// #include <RcppArmadillo.h>
using namespace Rcpp;

//' @title Calculate Multitaper Adaptive Weights
//' 
//' @description Used to calculate the adaptive weights on eigenspectra using a bound on the 
//' broadband bias.
//' 
//' @param eigenSpec A \code{matrix} with columns containing the eigenspectra
//' @param ev A \code{vector} with the eigenvalues of the Slepian tapers.
//' @param k A \code{numeric} integer representing the number of tapers to use.
//' @param nFFT A \code{numeric} integer representing the total number of frequency bins 
//' that were used to estimate the spectra.
//' 
//' @details When estimating the variance, this function will overestimate because 
//' it includes the zero and nyquist frequency twice (should only be once).
//' 
//' @return A \code{matrix} containing the weights for each eigenspectra at each 
//' frequency.
// [[Rcpp::export]]
NumericMatrix adaptiveWeightsCpp(NumericMatrix eigenSpec, NumericVector ev, int k, int nFFT) {
  int nfreq = nFFT/2 + 1;
  NumericVector sHat2(nfreq); // nfreq = nFFT/2 + 1
  NumericVector sHat(nfreq);
  NumericVector B(k); // broadband bias
  double sigmasq = 0.0; // variance
  double ssd; // sum of squared weights
  NumericMatrix d(nfreq, k); // weight matrix
  
  // pilot spectrum and initial variance estimate
  for (int i = 0; i < nfreq; i++){
    sHat2[i] = 0.5 * (eigenSpec(i, 0) + eigenSpec(i, 1));
    sigmasq += 2.0*sHat2[i];
  }
  
  sigmasq = sigmasq / nFFT; // normalize the variance properly
  
  sHat = 2.0 * sHat2; // Rcpp sugar... vectorization
  
  while(mean(abs(sHat2 - sHat) / sHat) > 0.0001){
    for (int m = 0; m < nfreq; m++){
      sHat[m] = sHat2[m];
    }
    
    // recalculate broadband bias (max)
    for (int i = 0; i < k; i++){
      B(i) = (1.0 - ev[i]) * sigmasq;
    }
    
    for (int i = 0; i < nfreq; i++){
      sHat2[i] = 0.0;
      ssd = 0.0;
      for (int j = 0; j < k; j++){
        d(i, j) = sqrt(ev[j]) * sHat[i] / (ev[j] * sHat[i] + B[j]); // calculate weight at freq i, taper j
        ssd += d(i, j); // sum of square weights (d)
        sHat2[i] += pow(std::abs(d(i, j)), 2.0) * eigenSpec(i, j); // numerator of adaptive weight formula (5.3 in '82)
      }
      sHat2[i] = sHat2[i] / ssd;
    }
   
   sigmasq = 2.0 * sum(sHat2) / nFFT;
 }
  
  return d;
}

// [[Rcpp::export]]
double absCplx(Rcomplex x){
  double y = sqrt(pow(x.r, 2.0) + pow(x.i, 2.0));
  return(y);
}

// [[Rcpp::export]]
ComplexMatrix coherencyOffsetCpp(ComplexMatrix ykx, ComplexMatrix yky
                                   , ComplexMatrix dx, ComplexMatrix dy
                                   , int nTaper, int nFreqRange, int nOffset
                                   , IntegerVector range2Start, IntegerVector freqRangeIdx
                                   , bool forwardCoh) {
  ComplexMatrix coh(nOffset, nFreqRange);
  Rcomplex Sx, Sy, numerator, denom, ssdx, ssdy, ssdxy, zero, dxdy;
  zero.r = 0.0;
  zero.i = 0.0;

  std::fill(coh.begin(), coh.end(), zero); // fill coh with 0's (complex)
  
  for(int os = 0; os < nOffset; os++){ // 0 to 1180
    // Rcout << "Offset: " << os << "\n"; // for testing purposes
    for (int f = 0; f < nFreqRange; f++){ // 0 to 1179
      numerator = zero;
      denom = zero;
      ssdx = zero;
      ssdy = zero;
      ssdxy = zero;
      Sx = zero;
      Sy = zero;
      
      for (int k = 0; k < nTaper; k++){
        if (forwardCoh){
          numerator = numerator + (dx(freqRangeIdx[f], k) * dy(range2Start[os]+f, k)) * (ykx(freqRangeIdx[f], k) * internal::complex__Conj(yky(range2Start[os]+f, k)));
        }
        else
        {
          numerator = numerator + (dx(freqRangeIdx[f], k) * dy(range2Start[os]+f, k)) * (ykx(freqRangeIdx[f], k) * yky(range2Start[os]+f, k));
        }
        ssdxy = ssdxy + dx(freqRangeIdx[f], k) * dy(range2Start[os]+f, k);
        Sx.r += pow(absCplx(dx(freqRangeIdx[f], k) * ykx(freqRangeIdx[f], k)), 2.0);
        Sy.r += pow(absCplx(dy(range2Start[os]+f, k)*yky(range2Start[os]+f, k)), 2.0);
        // Don't need this I don't think - replaced with Sx and Sy and denom below
        //denom.r += pow(absCplx(dx(freqRangeIdx[f], k) * ykx(freqRangeIdx[f], k)), 2.0) * pow(absCplx(dy(range2Start[os]+f, k)*yky(range2Start[os]+f, k)), 2.0);
        // Don't need these below ... I don't think, see DJT '82 - p. 1089
        // ssdx.r += pow(dx(freqRangeIdx[f], k).r, 2.0);
        // ssdy.r += pow(dy(range2Start[os]+f, k).r, 2.0);
      }
      // denom.r = sqrt( (Sx.r / ssdx.r) * (Sy.r / ssdy.r) );
      denom.r = sqrt( (Sx.r) * (Sy.r) );
      // don't need the normalizing factors? See DJT '82 - p. 1089
      // coh(os, f) = (numerator / ssdxy) / denom;
      coh(os, f) = (numerator) / denom;
    }
  }

  return(coh);
}