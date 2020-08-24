#include <RcppArmadillo.h>
#include <set>
#include "samplers.h"
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]



struct eKFOUT{
  double att; // a_t|t
  double at;  // a_t+1 (but this is done at the beginning of the KF
                 // instead of the end)
  double Ptt; // same deal
  double Pt;
  double lik;
  double pred; // f(at) to estimate yt
};

// [[Rcpp::export]]
eKFOUT ekfmm1step(
    arma::vec llt, int state, int t, arma::vec yyt,
    double sig2eps, arma::vec mus, arma::vec sig2eta,
    double a0, double P0){
  
  // Creates parameter matrices and does 1 step ekf
  // Specific to the simplified multiplicative model
  // 
  // States are: (1,1) (1,2) (1,4) (2,2) (2,3) (3,1) (3,3) (4,1), 
  //            (1,3), (2,1), (3,2) (same as with usual model)
  int nstates = 11;
  
  double ct;
  double Zt;
  double HHt = 0.0;
  double dt = 0.0;
  double Tt = 1.0;
  double GGt = sig2eps;
  
  double lt = llt(t);
  double yt = yyt(t);
  
  std::set<int> nonlin = {1, 3, 4, 6, 8, 10};
  std::set<int> acc = {4, 6, 8};
  
  if(state == 5 || state == 9){
    dt = mus(0);
    Tt = 0.0;
    HHt = sig2eta(0);
  }
  if(nonlin.find(state) != nonlin.end()){
    dt = log(lt);
    dt += acc.find(state) == acc.end() ? log(1+mus(1)) : log(1-mus(1));
    HHt = sig2eta(1);
  }
  
  // KF
  // Update a_t-1|t-1 to a_t using CURRENT state
  double atemp = dt + Tt * a0;
  double Ptemp = HHt + Tt * P0 * Tt;
  // Evaluate Taylor approximation of exp(x)
  Zt = exp(atemp);
  ct = Zt * (1.0 - atemp);
  if(state == 2){
    ct += mus(2);
    GGt += sig2eta(2);
  }
  // Make predictions using a_t
  double pred = ct + Zt * atemp;
  double vt = yt - pred;
  double Ft = GGt + Zt * Ptemp * Zt;
  double Ftinv = 1/Ft;
  double Kt = Ptemp * Zt * Ftinv;
  // Incorporate current information into a_t|t
  double a1 = atemp + Kt * vt;
  double P1 = Ptemp - Ptemp * Zt * Kt;
  
  // Calculate likelihood
  double ftdet = std::abs(Ftinv);
  ftdet = 1 / std::abs(ftdet);
  double mahalanobis = vt * Ftinv * vt;
  mahalanobis += log(2*M_PI) + log(ftdet);
  double lik = exp(-1.0*mahalanobis/2);
  eKFOUT output = {a1, atemp, P1, Ptemp, lik, pred};
  return output;
}








List edpf(arma::uvec currentStates, arma::colvec w, int N,
         arma::mat transProbs, arma::vec lt, int step,
         arma::vec a0, arma::vec P0,
         double sig2eps, arma::vec mus, arma::vec sig2eta, arma::vec yt){
  int npart = currentStates.size();
  int K = transProbs.n_cols;  //number of possible discrete states

  arma::mat a1(K, npart);
  arma::mat P1(K, npart);                   
  arma::mat lik(npart, K);
  arma::vec a11(K);
  arma::vec P11(K);
  for(int part=0; part < npart; part++){
    for(int k=0;  k < K; k++){   //both loops together have effect of for each particle
      eKFOUT ekfout = ekfmm1step(lt, k, step, yt, sig2eps, mus, sig2eta,
                                 a0(part), P0(part));
      double atmp = ekfout.att;
      double Ptmp = ekfout.Ptt;
      a11.col(k) = atmp;
      P11.col(k) = Ptmp;
      lik(part,k) = ekfout.lik;
    }
    a1.col(part) = a11;
    P1.col(part) = P11;
  }
  arma::mat testLik = lik;
  lik %= transProbs.rows(currentStates);
  lik.each_col() %= w;
  w = arma::vectorise(lik);
  w = arma::normalise(w,1);
  w = resampleSubOptimal(w, N);
  if( !arma::any(w)) return List::create(Named("BadPars") = 1);
  arma::uvec positive = arma::find(w);
  int nkeep = positive.size();
  arma::colvec newW = w(positive);
  arma::vec aout(nkeep);
  arma::vec Pout(nkeep);
  arma::uvec newstates(nkeep);
  arma::uvec oldstates(nkeep);
  int i=0;
  for(int k=0; k<K; k++){
    for(int part=0; part < npart; part++){
      if(w(part + k*npart)>0){
        newstates(i) = k;
        oldstates(i) = part;
        aout(i) = a1(k,part);
        Pout(i) = P1(k,part);
        i++;
      }
    }
  }
  return List::create(Named("BadPars") = 0,
                      Named("a1") = aout,
                      Named("P1") = Pout,
                      Named("oldstates") = oldstates,
                      Named("newstates") = newstates,
                      Named("newW") = newW,
                      Named("testLik") = testLik);
}


// [[Rcpp::export]]
double egetloglike(arma::vec lt, arma::uvec path, arma::mat y,
                   double sig2eps, arma::vec mus, arma::vec sig2eta,
                   arma::vec a0, arma::vec P0){
  
  arma::uword n = y.n_cols;
  
  double a0_ = a0(path(0));
  double P0_ = P0(path(0));
  
  double llik = 0;
  double liktmp;
  int s;
  
  for(int iter=0; iter<n; iter++){
    s = path(iter);
    eKFOUT step = ekfmm1step(lt, s, iter, y, sig2eps, mus, sig2eta,
                             a0_, P0_);
    a0 = step.att;
    P0 = step.Ptt;
    liktmp = step.lik;
    llik -= log(liktmp);
  }
  return llik;
}
    

 

List einitializeParticles(arma::vec w0, int N,  arma::vec lt, arma::vec yt,
                          double sig2eps, arma::vec mus, arma::vec sig2eta,
                          arma::vec a0, arma::vec P0){
  arma::uvec particles = arma::find(w0);
  int npart = particles.n_elem;
  if(npart == 0){return List::create(Named("BadPars") = 1);}
  arma::vec weights = w0(particles);
  arma::vec a11(npart);
  arma::vec P11(npart);
  arma::vec lik(npart);
  for(int i = 0; i < npart; i++){
    int p = particles(i);
    eKFOUT kfout = ekfmm1step(lt, p, 0, yt, sig2eps, mus, sig2eta,
                              a0(p), P0(p));
    a11(i) = kfout.att;
    P11(i) = kfout.Ptt;
    lik(i) = kfout.lik;
  }
  lik %= weights;
  arma::vec w = arma::normalise(lik, 1);
  w = resampleSubOptimal(w, N);
  if( !arma::any(w)) return List::create(Named("BadPars") = 1);
  arma::uvec nz = arma::find(w);
  arma::vec newW = w(nz);
  arma::uvec newstates = particles(nz);
  arma::mat aout = a11(nz);
  arma::mat Pout = P11(nz);
  return List::create(Named("BadPars") = 0,
                      Named("a1") = aout,
                      Named("P1") = Pout,
                      Named("newstates") = newstates,
                      Named("newW") = newW);
}


// [[Rcpp::export]]
List ebeamSearch(arma::vec lt, arma::vec w0,
                 double sig2eps, arma::vec mus, arma::vec sig2eta,
                 arma::vec a0, arma::vec P0, arma::mat y,
                 arma::mat transProbs, int N){
  
  arma::vec yt = y.row(0);
  arma::uword n = y.n_cols;
  arma::uword K = transProbs.n_cols;
  double maxpart = pow(K,n);
  if(maxpart > N) maxpart = N;
  a0.resize(maxpart);
  P0.resize(maxpart);
  
  arma::uvec particles(maxpart, arma::fill::zeros);
  arma::uvec filler = arma::regspace<arma::uvec>(0,K-1);
  arma::umat paths(maxpart, n, arma::fill::zeros);
  arma::colvec weights(maxpart,arma::fill::zeros);
  List step = einitializeParticles(w0, N, lt, yt, sig2eps, mus, sig2eta, a0, P0);
  
  int BP = step["BadPars"];
  if(BP) return List::create(Named("LastStep") = 1);
  arma::vec newW = step["newW"];
  arma::uword CurrentPartNum = newW.n_elem;
  weights.head(CurrentPartNum) = newW;
  arma::mat a1tmp = step["a1"];
  a0.head(CurrentPartNum) = a1tmp;
  arma::mat P1tmp = step["P1"];
  P0.head(CurrentPartNum) = P1tmp;
  arma::uvec newS = step["newstates"];
  particles.head(CurrentPartNum) = newS;
  paths(arma::span(0, CurrentPartNum - 1), 0) = newS;
  arma::uword iter = 1;
  while(iter < n){
    step = edpf(particles.head(CurrentPartNum), weights.head(CurrentPartNum), 
                maxpart, transProbs, lt, iter,
                a0.head(CurrentPartNum), P0.head_cols(CurrentPartNum),
                sig2, mus, sig2eta, yt);
    int BP = step["BadPars"];
    if(BP) break;
    arma::vec newW = step["newW"];
    CurrentPartNum = newW.n_elem;
    weights.head(CurrentPartNum) = newW;
    if(CurrentPartNum < maxpart) weights.tail(maxpart-CurrentPartNum) *= 0; // fix decreasing CurrentPartNum
    arma::mat a1tmp = step["a1"];
    a0.head(CurrentPartNum) = a1tmp;
    arma::mat P1tmp = step["P1"];
    P0.head(CurrentPartNum) = P1tmp;
    arma::uvec newS = step["newstates"];
    particles.head(CurrentPartNum) = newS;
    arma::uvec old = step["oldstates"];
    arma::umat tmp = paths.rows(old);
    paths.head_rows(CurrentPartNum) = tmp;   //could create issues if CurrentPartNum decreases over time?
    paths(arma::span(0,CurrentPartNum-1), iter) = newS;
    iter++;
  }
  return List::create(Named("paths") = paths,
                      Named("weights") = weights,
                      Named("LastStep") = iter++);
}

arma::mat ecreateTransMat(arma::vec transprobs){
  
  int nstates = 11;
  arma::mat transMat(nstates, nstates, arma::fill::zeros);
  if(transprobs.n_elem<6) transprobs.resize(6);
  transMat(0,0) = transprobs(0); // (1,1) -> (1,1)
  transMat(0,1) = transprobs(1); // (1,1) -> (1,2)
  transMat(0,2) = 1-transprobs(0)-transprobs(1)-transprobs(4); // (1,1) -> (1,4)
  transMat(0,8) = transprobs(4); // (1,1) -> (1,3)
  transMat(1,3) = 1; // (1,2) -> (2,2)
  transMat(2,7) = 1; // (1,4) -> (4,1)
  transMat(3,3) = transprobs(2); // (2,2) -> (2,2)
  transMat(3,4) = 1 - transprobs(2) - transprobs(5); // (2,2) -> (2,3)
  transMat(3,9) = transprobs(5); // (2,2) -> (2,1)
  transMat(4,6) = 1; // (2,3) -> (3,3)
  transMat(5,0) = 1; // (3,1) -> (1,1)
  transMat(6,5) = transprobs(3); // (3,3) -> (3,1)
  transMat(6,10) = transprobs(6); // (3,3) -> (3,2)
  transMat(6,6) = 1 - transprobs(3) - transprobs(6); // (3,3) -> (3,3)
  transMat(7,0) = 1; // (4,1) -> (1,1)
  transMat(8,6) = 1; // (1,3) -> (3,3)
  transMat(9,0) = 1; // (2,1) -> (1,1)
  transMat(10,3) = 1; // (3,2) -> (2,2)
  
  return transMat;
}





// [[Rcpp::export]]
List ekalman(arma::vec lt, arma::uvec path, arma::mat y,
             double sig2eps, arma::vec mus, arma::vec sig2eta,
             double a0, double P0){
  // What if I want different initial state (instead of 0)?
  

  arma::uword n = y.n_cols;
  
  // output storage
  arma::vec llik = arma::zeros(n);
  arma::vec at(n, arma::fill::zeros); // predicted mean E[a_t | y_1:t-1]
  arma::vec att(n, arma::fill::zeros); // updated mean E[a_t | y_1:t]
  arma::vec Pt(n,arma::fill::zeros); // predicted variance
  arma::vec Ptt(n,arma::fill::zeros); // predicted variance
  arma::vec preds(n,arma::fill::zeros); // predictions  
  arma::vec ahat(n,arma::fill::zeros); // smoothed mean E[a_t | y_1:n]
  arma::vec ests(n,arma::fill::zeros); //smoothed estimates
  double liktmp = 0.0;
  
  
  // Kalman filtering
  arma::uword iter=0;
  double a0_ = a0;
  double P0_ = P0;
  int s;
  while(iter < n){
    s = path(iter);
    eKFOUT step = ekfmm1step(lt, s, iter, y, sig2eps, mus, sig2eta, a0_, P0_);
  
    at(iter) = step.at;
    att(iter) = step.att;
    Pt(iter) = step.Pt;
    Ptt(iter) = step.Ptt;
    preds(iter) = step.pred;
    a0_ = step.att;
    P0_ = step.Ptt;
    liktmp = step.lik;
    llik(iter) += log(liktmp);
    iter++;
  }
  
  // Kalman smoothing. Based on Durbin and Koopman (4.70)
  iter--;
  ahat(iter) = att(iter);
  double Tt;
  double P00;
  double a00;
  
  while(iter > 0){
    s = path(iter);
    
    ests(iter) = exp(ahat(iter));
    if(s == 2) ests(iter) += mus(2);
    P00 = 1/Pt(iter);
    a00 = ahat(iter) - at(iter);
    
    iter--; // This is important, need T after the increment rather than before.
    s = path(iter);
    Tt = (s==5 || s==9) ? 0.0 : 1.0;
    ahat(iter) = att(iter) + Ptt(iter) * Tt * P00 * a00;
  }
  
  ests(0) = exp(ahat(0));
  if(s==2) ests(0) += mus(2);
    
  return List::create(Named("at") = at, Named("Pt") = Pt, Named("preds") = preds,
                            Named("ahat") = ahat,
                            Named("att") = att, Named("Ptt") = Ptt,
                            Named("ests") = ests, Named("llik") = llik);
}

