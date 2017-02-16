#include <RcppArmadillo.h>
#include <math.h>
#include <complex>
using namespace Rcpp;

const double pi = M_PI;

// Transfer function code - quick quick! (maybe.. )

//' Slow Fourier transform
//' 
//' Computes the Fourier transfrom of the data \code{x} at the frequency \code{f}.
//' 
//' @param x A complex vector to be Fourier transformed
//' @param f A \code{numeric} value indicating the frequency
//'  at which to calculate the Fourier tranform
//' @param time A vector containing the times (sampling rate sort of).
//' 
//' @details My understanding is that performing the Fourier transform in this way 
//' is numerically less stable than using an FFT and associated algorithms.
//' This is something that should be addressed by the author.
//' 
//' @return complex value of the Fourier transform of x at f.
//' 
//' @export
// [[Rcpp::export]]
arma::cx_double sftCpp(arma::cx_vec x, arma::cx_double f, arma::cx_vec time) {
  arma::cx_double ft(0.0, 0.0), eye(0.0, 1.0);
  int n = x.size();
  
  arma::cx_double negTwoPi(-2.0*pi, 0.0);
  
  for (int i = 0; i < n; i++){
    ft += (x[i] * std::exp(negTwoPi * eye * f * time[i]));
  }
  
  return ft;
}

//' Ordinary least squares frequency-domain regression
//' 
//' Estimates coefficients by frequency (transfer function) between 
//' time-domain inputs \code{x} and response \code{y}.
//' 
//' @param y a \code{list} of complex matrices (containing real data probably) 
//' - each matrix is a block in time and contains response * dpss - n rows, k columns
//' @param x is a \code{list} of \code{lists}: each sublist is a block of time containing 
//' numeric matrices like y - one matrix for each predictor
//' @param time A vector containing times at which the data are sampled for each block.
//' (you can just supply 1:n)
//' @param n An \code{integer} indicating the number of samples in each block.
//' @param npredictor An \code{integer} indicating the number of inputs (X's)
//' @param ntaper An \code{integer} indicating the number of tapers that were used.
//' @param freq is a list of frequencies at which to estimate the transfer function 
//' - someone needs to take into account that there will be a gap at beginning and 
//' end of size max(fOffset).
//' @param fOffset a vector indicating which frequencies of predictors to use 
//' (should at least be 0)
//' 
//' @details I need to revisit this and think of a way to NOT use the slow Fourier transform.
//' Also, probably want to be using an SVD regression to get the transfer function.
//' 
//' @export
// [[Rcpp::export]]
arma::cx_mat olsTf(List x, List y, arma::cx_vec time, int n
                     , int npredictor, int ntaper
                     , arma::cx_vec freq, arma::cx_vec fOffset){
  
  int nblock = x.size(), nfreq = freq.size(), noffset = fOffset.size(), cur, zCur;
  arma::cx_mat Z(nfreq, npredictor), X(nblock*ntaper*noffset, npredictor * noffset);
  arma::cx_vec Y(nblock*ntaper*noffset), Ztmp(npredictor);
  
  // initialize these?
  arma::cx_mat Ytmp = as<arma::cx_mat>(y[1]);
  arma::cx_mat Xtmp = as<arma::cx_mat>(as<List>(x[1])[1]);
  
  if (nblock != y.size()) stop("Must have the same number of blocks.");
  
  // zCur = 0;
  for (int f = 0; f < nfreq; f++){
    
    // The next 4 for-loops set up the response vector and design matrix
    cur = 0;
    for (int b = 0; b < nblock; b++){
      // Rcpp::ComplexMatrix Ytmp = as<Rcpp::ComplexMatrix>(y[b]);
      Ytmp = as<arma::cx_mat>(y[b]);
      
      for (int k = 0; k < ntaper; k++){
        Y(cur) = sftCpp(Ytmp.col(k), freq(f), time);
        for (int s = 0; s < noffset; s++){
          if (s > 0){
            Y(cur) = Y(cur-1); // respones variable keeps the same value
          }
          
          for (int p = 0; p < npredictor; p++){
            Xtmp = as<arma::cx_mat>(as<List>(x[b])[p]);
            X(cur, p) = sftCpp(Xtmp.col(k), freq(f) + fOffset(s), time);
          }
          cur++;
        }
      }
    }
    
    // do the actual least squares beta = (t(Conj(X)) * X)^(-1) * t(Conj(X)) * Y
    Ztmp = (arma::inv(arma::trans(X) * X) * arma::trans(X)) * Y;
    
    for (int i = 0; i < Ztmp.n_rows; i++){
      // Z(zCur, i) = Ztmp(i);
      Z(f, i) = Ztmp(i);
    }
    // zCur++;
  }
  
  return Z;
}


//' Ordinary least squares frequency-domain regression on eigencoefficients
//' 
//' Estimates coefficients by frequency (transfer function) between 
//' inputs \code{x} and response \code{y}.
//' 
//' @param x \code{list} of inputs that is a result of \code{taper()}.
//' @param y \code{list} of responses that is a result of \code{taper()}.
//' @param freqIdx \code{numeric} vector indicating at which frequencies the 
//' coefficients should be estimated
//' @param fOffsetIdxLst a \code{list} containing what offset frequencies should be used 
//' (in a standard transfer function, this would contain the zero frequency offset)
//' 
//' @details Note to me: I should look at exactly how this is working (re: fOffsetIdxLst)
//' 
//' @return a \code{matrix} with columns the transfer function relating to the input 
//' that was in the same column.
//' 
//' @export
// [[Rcpp::export]]
arma::cx_mat olsTfEigen(List x, List y
                          , arma::ivec freqIdx, List fOffsetIdxLst){
  
  int nblock = x.size(), nfreq = freqIdx.size(), totalnoffset = 0, cur, curPred;
  int npredictor = as<List>(x[1]).size();
  arma::ivec noffset = arma::zeros<arma::ivec>(fOffsetIdxLst.size());
  for (int i = 0; i < fOffsetIdxLst.size(); i++){
    noffset(i) = as<arma::ivec>(fOffsetIdxLst[i]).size();
    totalnoffset += noffset(i);
  }
  
  // initialize these? - YES
  arma::cx_mat Ytmp = as<arma::cx_mat>(y[1]);
  int ntaper = Ytmp.n_cols;
  arma::cx_mat Xtmp = as<arma::cx_mat>(as<List>(x[1])[1]);
  arma::cx_mat Z(nfreq, totalnoffset), X(nblock*ntaper, totalnoffset);
  arma::cx_vec Y(nblock*ntaper), Ztmp(totalnoffset);
  
  if (nblock != y.size()) stop("Must have the same number of blocks.");
  
  // zCur = 0;
  for (int f = 0; f < nfreq; f++){
    // The next 4 for-loops set up the response vector and design matrix
    cur = 0;
    for (int b = 0; b < nblock; b++){
      // Rcpp::ComplexMatrix Ytmp = as<Rcpp::ComplexMatrix>(y[b]);
      Ytmp = as<arma::cx_mat>(y[b]);
      
      for (int k = 0; k < ntaper; k++){
        Y(cur) = Ytmp(freqIdx(f), k);
        curPred = 0;
        for (int p = 0; p < as<List>(x[1]).size(); p++){
          Xtmp = as<arma::cx_mat>(as<List>(x[b])[p]);
          arma::ivec fOffsetIdx = as<arma::ivec>(fOffsetIdxLst[p]);
          for (int s = 0; s < fOffsetIdx.size(); s++){
            X(cur, curPred) = Xtmp(freqIdx(f) + fOffsetIdx(s), k);
            curPred++;
          }
        }
        cur++;
      }
    }
    
    // do the actual least squares beta = (t(Conj(X)) * X)^(-1) * t(Conj(X)) * Y
    Ztmp = (arma::inv(arma::trans(X) * X) * arma::trans(X)) * Y;
    
    for (int i = 0; i < Ztmp.n_rows; i++){
      // Z(zCur, i) = Ztmp(i);
      Z(f, i) = Ztmp(i);
    }
    // zCur++;
  }
  
  return Z;
}