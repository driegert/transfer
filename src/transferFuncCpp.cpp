// #include <Rcpp.h> // not needed I guess...
#include <RcppArmadillo.h>
#include <math.h>
#include <complex>
using namespace Rcpp;

const double pi = M_PI;

// Transfer function code - quick quick! (maybe.. )

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

/* 
 * y is a list of complex matrices (containing real data probably)
 * - each matrix is a block in time and contains 
 * response * dpss - n rows, k columns
 */
/*
 *  x is a list of lists:
 *  each sublist is a block of time containing
 *  numeric matrices like y - one matrix for each predictor
 */
/*
 *  freq is a list of frequencies at which to estimate the transfer function
 *  - someone needs to take into account that there will be a gap at beginning and 
 *  end of size max(fOffset).
 */
// fOffset tells us which frequencies of predictors to use (should at least be 0)
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



// BACKUP FOR NOW - JULY 29 @ 9:42am
// ComplexMatrix olsTf(List x, List y, int n
//                       , int npredictor, int ntaper
//                       , ComplexVector freq, ComplexVector fOffset){
//   
//   int nblock = x.size(), nfreq = freq.size(), noffset = fOffset.size(), cur;
//   ComplexMatrix Z(nfreq, npredictor), X(nblock*ntaper*noffset, npredictor);
//   ComplexVector Y(nblock*ntaper*noffset);
//   
//   if (nblock != y.size()) stop("Must have the same number of blocks.");
//   
//   for (int f = 0; f < nfreq; f++){
//     
//     // The next 4 for-loops set up the response vector and design matrix
//     cur = 0;
//     for (int b = 0; b < nblock; b++){
//       Rcpp::ComplexMatrix Ytmp = as<Rcpp::ComplexMatrix>(y[b]);
//       
//       for (int k = 0; k < ntaper; k++){
//         Y(cur) = sftCpp(Ytmp(_, k), freq(f));
//         for (int s = 0; s < noffset; s++){
//           Y(std::min(0, cur)) = Y(cur-1); // respones variable keeps the same value
//           for (int p = 0; p < npredictor; p++){
//             Rcpp::ComplexMatrix Xtmp = as<Rcpp::ComplexMatrix>(as<List>(x[b])[p]);
//             X(cur, p) = sftCpp(Xtmp(_, k), freq(f) + fOffset(s));
//           }
//           cur++;
//         }
//       }
//     }
//     
//     // do the actual least squares beta = (t(Conj(X)) * X)^(-1) * t(Conj(X)) * Y
//     
//     arma::cx_mat XX = as<arma::cx_mat>(X);
//     arma::cx_mat YY = as<arma::cx_mat>(Y);
//     
//     arma::cx_mat ZZ = (arma::inv((arma::trans(XX)) * XX) * arma::trans(XX)) * YY;
//     
//     for (int j = 0; j < npredictor; j++){
//       Z(f, j) = as<Rcomplex>(ZZ[j, 1]);
//     }
//   }
//   
//   return Z;
// }


// Rcomplex sftCpp(ComplexVector x, Rcomplex f) {
//   Rcomplex expFt;
//   // ComplexVector ft(f.size());
//   Rcomplex ft;
//   ft.r = 0;
//   ft.i = 0;
//   int n = x.size();
//   
//   // for (int i = 0; i < ft.size(); i++){
//   //   ft(i).r = 0;
//   //   ft(i).i = 0;
//   // }
//   
//   // for (int j = 0; j < ft.size(); j++){
//   for (int i = 0; i < n; i++){
//     expFt.r = cos(2 * pi * f.r * i);
//     expFt.i = -sin(2 * pi * f.r * i);
//     ft = ft + (x(i) * expFt);
//   }
//   //}
//   
//   return ft;
// }