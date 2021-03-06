#ifndef __SAMPLERS_H
#define __SAMPLERS_H


#include <RcppArmadillo.h>
using namespace Rcpp;


arma::uvec SampleNoReplace(arma::uvec x, int size) {
  int nOrig = x.size();
  arma::uvec index(size);
  arma::uvec sub(nOrig);
  for (int ii = 0; ii < nOrig; ii++) sub(ii) = ii;
  RNGScope scope4;
  for (int ii = 0; ii < size; ii++) {
    int jj = floor(nOrig * runif(1)[0]);
    index(ii) = sub(jj);
    // replace sampled element with last, decrement
    sub(jj) = sub(--nOrig);
  }
  arma::uvec ret(size);
  ret = x(index);
  return(ret);
}


arma::vec resampleSubOptimal(arma::vec w, int N){
  int M = w.size();
  double tol = 1e-10;
  arma::vec ws = w;
  ws.elem( arma::find( ws < tol) ).zeros();
  arma::vec nzz = arma::nonzeros(ws);
  arma::uword nz = nzz.n_elem;
  if(M <= N || nz <= N){
    return ws;
  }
  
  arma::vec w1 = ws;
  typedef std::vector<double> stdvec;
  stdvec z = arma::conv_to<stdvec>::from(w1);
  std::nth_element(z.begin(), z.end()-N, z.end());
  arma::vec z1 = arma::conv_to<arma::vec>::from(z);
  double minprob = z1(M-N);                              //This is the Nth largest element of z1.
  
  ws.elem(arma::find(ws < minprob)).zeros();
  arma::uvec ties = arma::find(ws==minprob);
  arma::uvec keep = arma::find(ws > minprob);
  int nkeep = keep.size();
  int tt = ties.size();
  if(tt > 1){
    arma::vec probs(tt);
    double p = 1.0 / tt;
    probs.fill(p);                                    //MICHAEL: probs is never used?
    int dontneed = tt - (N - nkeep);
    arma::uvec idx = SampleNoReplace(ties, dontneed);
    ws.elem(idx).zeros();
  }
  ws = arma::normalise(ws,1);
  return ws;
}


arma::colvec resampleOptimal(arma::colvec w, int N){
  // no zeros no dups?? unused, doesn't seem to work
  int M = w.size();
  double tol = 1e-10;
  arma::colvec ws = w;
  ws.elem( arma::find( ws < tol) ).zeros();
  arma::vec nzz = arma::nonzeros(ws);
  arma::uword nz = nzz.n_elem;
  if(M <= N || nz <= N){ 
    return ws;
  }
  
  // Hard case.
  ws = arma::sort(ws);
  int Ak = 0;
  double Bk = 1.0;
  int i = M;
  while(i > M - N){
    i--;
    if(ws(i) < tol || Bk <= tol){
      w.elem( arma::find( w < tol) ) .zeros();
      return w;
    }
    if(Bk/ws(i) + Ak >= N) break;
    Ak++;
    Bk -= ws(i);
  }
  double cinv = Bk / (N-Ak);
  
  // Set 1
  arma::vec NewW;
  NewW.zeros(M);
  int L=0;
  double K=0;
  for(i=0; i<M; i++){
    if(cinv < w(i)){
      NewW(i) += w(i);
      L++;
    }else{
      K += w(i);  
    }
  }
  K /= (N-L);
  // Set 2
  RNGScope scope;
  double U1 = runif(1)[0] * K;
  for(int i=0; i<M; i++){
    if(NewW(i)==0){
      U1 -= w(i);
      if(U1 < 0){
        NewW(i) += cinv;
        U1 += K;
      }
    }
  }
  return NewW;
}

#endif 