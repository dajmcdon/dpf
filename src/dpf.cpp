#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]


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

//[[Rcpp::export]]
arma::vec resampleSubOptimal(arma::vec w, int N){
  int M = w.size();
  double tol = 1e-10;
  arma::vec ws = w;
  ws.elem( arma::find( ws < tol) ).zeros();
  // ws = arma::normalise(ws,1);
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

//[[Rcpp::export]]
arma::colvec resampleOptimal(arma::colvec w, int N){
  // no zeros no dups?? unused, doesn't seem to work
  int M = w.size();
  double tol = 1e-10;
  arma::colvec ws = w;
  ws.elem( arma::find( ws < tol) ).zeros();
  // ws = arma::normalise(ws,1);
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

struct KFOUT{
  arma::mat att; // a_t|t
  arma::mat at;  // a_t+1 (but this is done at the beginning of the KF
                 // instead of the end)
  arma::mat Ptt; // same deal
  arma::mat Pt;
  double lik;
  arma::mat pred; // f(at) to estimate yt
};

    

KFOUT kf1step(arma::mat a0, arma::mat P0, arma::mat dt,
             arma::mat ct, arma::mat Tt,
             arma::mat Zt, arma::mat HHt, arma::mat GGt, arma::mat yt) {
  // Reshape matrices
  int m = a0.size();//number of elements.
  int d = yt.size();
  P0.reshape(m,m);
  Tt.reshape(m,m);
  Zt.reshape(d,m);
  HHt.reshape(m,m); // State var (see fkf vs Durbin)
  GGt.reshape(d,d); // Obs var
  arma::mat Ptemp(m,m);
  arma::mat atemp(m,1);

  // KF
  // Update a_t-1|t-1 to a_t using CURRENT state
  atemp = dt + Tt * a0;
  Ptemp = HHt + Tt * P0 * Tt.t();
  
  // Make predictions using a_t
  arma::mat pred = ct + Zt * atemp;
  arma::mat vt = yt - pred;
  arma::mat Ft = GGt + Zt * Ptemp * Zt.t();
  arma::mat Ftinv = arma::inv(Ft);
  arma::mat Kt = Ptemp * Zt.t() * Ftinv;
  
  // Incorporate current information into a_t|t
  arma::mat a1 = atemp + Kt * vt;
  arma::mat P1 = Ptemp - Ptemp * Zt.t() * Kt.t();
  // a1 = dt + Tt * a1;
  // P1 = HHt + Tt * P1 * Tt.t();

  // Calculate likelihood
  double ftdet = arma::det(Ftinv);
  ftdet = 1 / std::abs(ftdet);
  double mahalanobis = arma::as_scalar(vt.t() * Ftinv * vt);
  mahalanobis += yt.size()* log(2*M_PI) + log(ftdet);
  double lik = exp(-1.0*mahalanobis/2);
  P1.reshape(m*m, 1);
  KFOUT output = {a1, atemp, P1, Ptemp, lik, pred};
  return output;
  /*return List::create(Named("a1") = a1,
                      Named("P1") = P1,
                      Named("lik") = lik,
                      Named("pred") = pred);*/
}

// [[Rcpp::export]]
List kf1stepR(arma::mat a0, arma::mat P0, arma::mat dt,
              arma::mat ct, arma::mat Tt,
              arma::mat Zt, arma::mat HHt, arma::mat GGt, arma::mat yt) {
    KFOUT output = kf1step(a0, P0, dt, ct, Tt, Zt, HHt, GGt, yt);
    
    return List::create(Named("att") = output.att,
                        Named("at") = output.at,
                        Named("Ptt") = output.Ptt,
                        Named("Pt") = output.Pt,
                        Named("lik") = output.lik,
                        Named("pred") = output.pred);
}


  

KFOUT ks1step(arma::mat r1, arma::mat N1, 
              arma::mat pred, arma::mat Pt, 
              arma::mat ct, arma::mat Tt,
              arma::mat Zt, arma::mat GGt, arma::mat yt) {
    
    // Reshape matrices
    int m = r1.size();//number of elements.
    int d = yt.size();
    Pt.reshape(m,m);
    Tt.reshape(m,m);
    Zt.reshape(d,m);
    GGt.reshape(d,d); // Obs var
    
    // KS
    // arma::mat pred = ct + Zt * at;
    arma::mat vt = yt - pred;
    arma::mat Ft = GGt + Zt * Pt * Zt.t();
    arma::mat Ftinv = arma::inv(Ft);
    arma::mat Lt = Tt - Tt * Pt * Zt.t() * Ftinv * Zt;
    arma::mat r0 = Zt.t() * Ftinv * vt + Lt.t() * r1;
    arma::mat N0 = Zt.t() * Ftinv * Zt + Lt.t() * N1 * Lt;
    
    KFOUT output = {r0, r0, N0, N0, 1, r0}; // Note: we don't need anything but r0 and N0
    return output;
}


//[[Rcpp::export]]
arma::mat HHcreate(arma::mat Rt, arma::mat Qt, int r, int q){
  arma::uword K = Rt.n_cols;
  arma::mat Rtmp(r,q);
  arma::mat Qtmp(q,q);
  arma::mat HHt(r*r, K);
  for(arma::uword iter=0; iter < K; iter++){
    Rtmp = reshape(Rt.col(iter), r, q);
    Qtmp = reshape(Qt.col(iter), q, q);
    HHt.col(iter) = arma::vectorise(Rtmp * Qtmp * Rtmp.t());
  }
  return(HHt);
}

//currentStates: vector of the current discrete state for each particle
//w: resampling weight for each particle
//N: max number of particles
//transProbs: matrix of transition probabilities for the discrete states
//a0: estimated means of the current continous state for each particle
//P0: estimated variances of the continous state for each particle
//dt: continuous state intercept. each column is an intercept
//ct: observation (yt) intercept
//Tt: continous state slope
//Zt: observation (yt) slope
//HHt: variance of the predicted mean of the continuous state
//GGt: variance of the observation
//yt: observation
List dpf(arma::uvec currentStates, arma::colvec w, int N,                      //MICHAEL: What are these?
         arma::mat transProbs,
         arma::mat a0, arma::mat P0,
         arma::mat dt, arma::mat ct, arma::mat Tt, arma::mat Zt, //Matrices stored in columns of these variables. Re turned into matrices in other functions.
         arma::mat HHt, arma::mat GGt, arma::vec yt){
  int npart = currentStates.size();
  int m = a0.n_rows;                            //dimension of each particle (the continuous part)
  int K = dt.n_cols;                            //number of possible discrete states
  int mm = m*m;
  arma::cube a1(m, K, npart);
  arma::cube P1(mm, K, npart);                   
  arma::mat lik(npart, K);
  arma::mat a11(m, K);
  arma::mat P11(mm, K);
  for(int part=0; part < npart; part++){
    for(int k=0;  k < K; k++){   //both loops together have effect of for each particle
      KFOUT kfout = kf1step(a0.col(part), P0.col(part), dt.col(k), ct.col(k),
                                       Tt.col(k), Zt.col(k), HHt.col(k), GGt.col(k),
                                       yt);
      arma::mat atmp = kfout.att;
      arma::mat Ptmp = kfout.Ptt;
      a11.col(k) = atmp;
      P11.col(k) = Ptmp;
      lik(part,k) = kfout.lik;
    }
    a1.slice(part) = a11;
    P1.slice(part) = P11;
  }
  // arma::mat blah = transProbs.rows(currentStates);
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
  arma::mat aout(m, nkeep);
  arma::mat Pout(mm, nkeep);
  arma::uvec newstates(nkeep);
  arma::uvec oldstates(nkeep);
  int i=0;
  for(int k=0; k<K; k++){
    for(int part=0; part < npart; part++){
      if(w(part + k*npart)>0){
        newstates(i) = k;
        oldstates(i) = part;
        arma::mat a1tmp = a1.subcube(0,k,part,m-1,k,part);
        arma::mat P1tmp = P1.subcube(0,k,part,mm-1,k,part);
        aout.col(i) = a1tmp;
        Pout.col(i) = P1tmp;
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


//' Evaluate the likelihood given parameters and discrete states
//' 
//' @param pmats a list of parameter matrices for the kalman filter
//' @param path the desired path for hidden discrete states
//' @param y observations, each time point in a column
//' 
//' @return the negative log-likelihood
//' 
//' @export 
// [[Rcpp::export]]
double getloglike(List pmats, arma::uvec path, arma::mat y){
  arma::mat a0 = pmats["a0"];
  arma::mat P0 = pmats["P0"];
  arma::cube dt = pmats["dt"];
  arma::cube ct = pmats["ct"];
  arma::cube Tt = pmats["Tt"];
  arma::cube Zt = pmats["Zt"];
  arma::cube Rt = pmats["Rt"];
  arma::cube Qt = pmats["Qt"];
  arma::cube GGt = pmats["GGt"];
  
  arma::uword n = y.n_cols;
  arma::uword m = a0.n_rows;
  arma::uword mm = m*m;
  arma::uword d = y.n_rows;
  arma::uword dm = d*m;
  arma::uword dd = d*d;
  
  arma::uword dtvar = dt.n_slices > 1;
  arma::uword ctvar = ct.n_slices > 1;
  arma::uword Ttvar = Tt.n_slices > 1;
  arma::uword Ztvar = Zt.n_slices > 1;
  arma::uword Rtvar = Rt.n_slices > 1;
  arma::uword Qtvar = Qt.n_slices > 1;
  arma::uword GGtvar = GGt.n_slices > 1;
  arma::mat HHt(m,m);
  arma::mat R(mm,1);
  arma::mat Q(mm,1);
  
  arma::mat aa0 = a0.col(path(0));
  arma::mat PP0 = reshape(P0.col(path(0)), m, m);
  
  double llik = 0;
  double liktmp;
  
  for(arma::uword iter=0; iter<n; iter++){
    arma::uword s = path(iter);
    if(iter==0 || Rtvar || Qtvar){ 
      R = Rt.subcube(0,s,iter*Rtvar,arma::size(mm,1,1));
      Q = Qt.subcube(0,s,iter*Qtvar,arma::size(mm,1,1));
      R.reshape(m,m);
      Q.reshape(m,m);
      HHt = R * Q * R.t();
    }
    KFOUT step = kf1step(aa0, PP0, dt.subcube(0,s,iter*dtvar,arma::size(m,1,1)),
                                    ct.subcube(0,s,iter*ctvar,arma::size(d,1,1)), 
                                    Tt.subcube(0,s,iter*Ttvar,arma::size(mm,1,1)),
                                    Zt.subcube(0,s,iter*Ztvar,arma::size(dm,1,1)), 
                                    HHt, GGt.subcube(0,s,iter*GGtvar,arma::size(dd,1,1)), 
                                    y.col(iter));
    
    aa0 = step.att;
    PP0 = step.Ptt;
    liktmp = step.lik;
    llik -= log(liktmp);
  }
  return llik;
}
    

 
//' Create parameter matrices as in Gu (2017)
//' 
//' @param lt durations between successive notes in the score
//' @param sig2eps variance of the observation noise
//' @param mus vector of 4 mean parameters
//' @param vector of 4 state variance parameters
//' @param vector of 4 transition probabilities
//' 
//' @return List with components as appropriate for Kalman filtering or Beam Search
//' 
//' @export    
// [[Rcpp::export]]
List yupengMats(arma::vec lt, double sig2eps, arma::vec mus,
                arma::vec sig2eta, arma::vec transprobs){ 
  //confirm that t's stay in same order for each matrix
  // in each section, we have:
  //   3 state means (tempo, accel, stress), 3 state variances (same), 1 obs variance
  //   4 transition matrix parameters
  // up to 22 parameters per performance
  // for now, everything constant, except mu_tempo
  int nstates = 8;
  int d = 1;
  int m = 2;
  int mm = m*m;
  int n = lt.n_elem;
  arma::mat a0(m, nstates, arma::fill::zeros);
  a0.row(0) += mus(0);
  a0(1,4) += mus(1);
  a0(1,6) += mus(1);
  arma::mat P0(mm, nstates, arma::fill::zeros);
  P0.row(0) += sig2eta(0);
  P0(3,4) += sig2eta(1);
  P0(3,6) += sig2eta(1);
  arma::cube dt(m, nstates, n, arma::fill::zeros);
  dt.tube(0,1) = lt*mus(1);
  dt.tube(0,4) = -lt*mus(1);
  dt.tube(1,1) += mus(1);
  dt.tube(1,4) -= mus(1);
  dt.tube(1,2) += mus(2);
  arma::vec tempo_mus(n, arma::fill::zeros);
  // tempo_mus.elem( find(temposwitch > 0) ) += mus(1);
  //tempo_mus.elem( find(temposwitch == 0) ) += mus(0);
  dt.tube(0,5) += mus(0);
  arma::cube ct(d, nstates, 1, arma::fill::zeros);
  arma::cube Tt(mm, nstates, n, arma::fill::zeros);
  Tt.tube(0,0,0,4) += 1;
  Tt.tube(0,6,0,7) += 1;
  Tt.tube(3,3) += 1;
  Tt.tube(3,6) += 1;
  Tt.tube(2,6) = lt;
  Tt.tube(2,3) = lt;
  arma::cube Zt(m, nstates, 1, arma::fill::zeros);
  Zt.tube(0,0,0,7) += 1;
  Zt(1,2,0) += 1;
  arma::cube Rt(mm, nstates, n, arma::fill::zeros);
  // Rt.tube(0,0) += 1;
  Rt.tube(0,0,0,7) += 1;
  // Rt.tube(0,5) += 1;
  Rt.tube(3,1,3,2) += 1;
  Rt.tube(3,4) += 1;
  Rt.tube(2,1) = lt;
  Rt.tube(2,4) = lt;
  arma::cube Qt(mm, nstates, 1, arma::fill::zeros);
  Qt.tube(0,0,0,7) += sig2eta(0);
  Qt.tube(3,0,3,7) += sig2eta(1);
  Qt.tube(3,2) -= sig2eta(1) - sig2eta(2);
  Qt.tube(0,5) = sig2eta(3);
  arma::cube GGt(d, nstates, 1, arma::fill::ones);
  GGt *= sig2eps;
  arma::mat transMat(nstates, nstates, arma::fill::zeros);
  transMat(0,0) = transprobs(0); // (1,1) -> (1,1)
  transMat(0,1) = transprobs(1); // (1,1) -> (1,2)
  transMat(0,2) = 1-transprobs(0)-transprobs(1); // (1,1) -> (1,4)
  transMat(1,3) = transprobs(2); // (1,2) -> (2,2)
  transMat(1,4) = 1-transprobs(2); // (1,2) -> (2,3)
  transMat(2,7) = 1; // (1,4) -> (4,1)
  transMat(3,3) = transprobs(2); // (2,2) -> (2,2)
  transMat(3,4) = 1-transprobs(2); // (2,2) -> (2,3)
  transMat(4,5) = transprobs(3); // (2,3) -> (3,1)
  transMat(4,6) = 1-transprobs(3); // (2,3) -> (3,3)
  transMat(5,0) = transprobs(0); // (3,1) -> (1,1)
  transMat(5,1) = transprobs(1); // (3,1) -> (1,2)
  transMat(5,2) = 1-transprobs(0)-transprobs(1); // (3,1) -> (1,4)
  transMat(6,5) = transprobs(3); // (3,3) -> (3,1)
  transMat(6,6) = 1-transprobs(3); // (3,3) -> (3,3)
  transMat(7,0) = transprobs(0); // (4,1) -> (1,1)
  transMat(7,1) = transprobs(1); // (4,1) -> (1,2)
  transMat(7,2) = 1-transprobs(0)-transprobs(1); // (4,1) -> (1,4)
  return List::create(Named("a0") = a0, Named("P0") = P0,
                      Named("dt") = dt, Named("ct") = ct,
                      Named("Tt") = Tt, Named("Zt") = Zt,
                      Named("Rt") = Rt, Named("Qt") = Qt,
                      Named("GGt") = GGt,
                      Named("transMat") = transMat);
}

List initializeParticles(arma::vec w0, int N, arma::mat a0, arma::mat P0,
                arma::mat dt, arma::mat ct, arma::mat Tt, arma::mat Zt,
                arma::mat HHt, arma::mat GGt, arma::vec yt){
    arma::uvec particles = arma::find(w0);
    int npart = particles.n_elem;
    if(npart == 0){return List::create(Named("BadPars") = 1);}
    arma::vec weights = w0(particles);
    arma::mat a11(a0.n_rows, npart);
    arma::mat P11(P0.n_rows, npart);
    arma::vec lik(npart);
    for(int i = 0; i < npart; i++){
        int p = particles(i);
        KFOUT kfout = kf1step(a0.col(p), P0.col(p), dt.col(p), ct.col(p),
                              Tt.col(p), Zt.col(p), HHt.col(p), GGt.col(p), yt);
        arma::mat atmp = kfout.att;
        arma::mat Ptmp = kfout.Ptt;
        a11.col(i) = atmp;
        P11.col(i) = Ptmp;
        lik(i) = kfout.lik;
    }
    lik %= weights;
    arma::vec w = arma::normalise(lik, 1);
    w = resampleSubOptimal(w, N);
    if( !arma::any(w)) return List::create(Named("BadPars") = 1);
    arma::uvec nz = arma::find(w);
    arma::vec newW = w(nz);
    arma::uvec newstates = particles(nz);
    arma::mat aout = a11.cols(nz);
    arma::mat Pout = P11.cols(nz);
    return List::create(Named("BadPars") = 0,
                        Named("a1") = aout,
                        Named("P1") = Pout,
                        Named("newstates") = newstates,
                        Named("newW") = newW);
}

//' Greedy HMM estimation given continuous hidden states
//' 
//' @param a0 a px1 matrix of state prior means
//' @param P0 a pxp matrix of state
//' @param dt dx1 or dxn matrix of 
//' 
//' @return List with components "paths", "weights", and "LastStep" 
//' 
//' @export 
// [[Rcpp::export]]
List beamSearch(arma::mat a0, arma::mat P0, arma::vec w0,
                arma::cube dt, arma::cube ct, arma::cube Tt, arma::cube Zt,
                arma::cube Rt, arma::cube Qt, arma::cube GGt, arma::mat yt,
                arma::mat transProbs, int N){
  arma::uword n = yt.n_cols;
  arma::uword K = transProbs.n_cols;
  double maxpart = pow(K,n);
  if(maxpart > N) maxpart = N;
  int m = a0.n_rows;
  a0.resize(m, maxpart);
  P0.resize(m*m, maxpart);
  // Determine if any parameter matrices are time varying. Assumes that each has either:
  // (1) 1 column per state or (2) n columns. Does not perform any checks.
  arma::uword dtvar = dt.n_slices > 1;
  arma::uword ctvar = ct.n_slices > 1;
  arma::uword Ttvar = Tt.n_slices > 1;
  arma::uword Ztvar = Zt.n_slices > 1;
  arma::uword Rtvar = Rt.n_slices > 1;
  arma::uword Qtvar = Qt.n_slices > 1;
  arma::uword GGtvar = GGt.n_slices > 1;
  arma::uword q = sqrtf(Qt.n_rows);
  arma::uvec particles(maxpart, arma::fill::zeros);
  arma::uvec filler = arma::regspace<arma::uvec>(0,K-1);
  arma::umat paths(maxpart, n, arma::fill::zeros);
  arma::colvec weights(maxpart,arma::fill::zeros);
  arma::mat HHt(m*m, K);
  HHt = HHcreate(Rt.slice(0), Qt.slice(0), m, q);
  List step = initializeParticles(w0, N, a0, P0, dt.slice(0), ct.slice(0),
                                  Tt.slice(0), Zt.slice(0), HHt, GGt.slice(0), yt.col(0));
  int BP = step["BadPars"];
  if(BP) return List::create(Named("LastStep") = 1);
  arma::vec newW = step["newW"];
  arma::uword CurrentPartNum = newW.n_elem;
  weights.head(CurrentPartNum) = newW;
  arma::mat a1tmp = step["a1"];
  a0.head_cols(CurrentPartNum) = a1tmp;
  arma::mat P1tmp = step["P1"];
  P0.head_cols(CurrentPartNum) = P1tmp;
  arma::uvec newS = step["newstates"];
  particles.head(CurrentPartNum) = newS;
  paths(arma::span(0, CurrentPartNum - 1), 0) = newS;
  arma::uword iter = 1;
  while(iter < n){
    if(Rtvar || Qtvar){ 
      HHt = HHcreate(Rt.slice(iter * Rtvar), Qt.slice(iter * Qtvar), m, q);
    }
    step = dpf(particles.head(CurrentPartNum), weights.head(CurrentPartNum), 
                              maxpart, transProbs, 
                              a0.head_cols(CurrentPartNum), P0.head_cols(CurrentPartNum),
                              dt.slice(iter * dtvar), ct.slice(iter * ctvar), 
                              Tt.slice(iter * Ttvar), Zt.slice(iter * Ztvar), 
                              HHt, GGt.slice(iter * GGtvar), yt.col(iter));
    int BP = step["BadPars"];
    if(BP) break;
    arma::vec newW = step["newW"];
    CurrentPartNum = newW.n_elem;
    weights.head(CurrentPartNum) = newW;
    arma::mat a1tmp = step["a1"];
    a0.head_cols(CurrentPartNum) = a1tmp;
    arma::mat P1tmp = step["P1"];
    P0.head_cols(CurrentPartNum) = P1tmp;
    arma::uvec newS = step["newstates"];
    particles.head(CurrentPartNum) = newS;
    arma::uvec old = step["oldstates"];
    arma::umat tmp = paths.rows(old);
    paths.head_rows(CurrentPartNum) = tmp;   //could create issues if CurrentPartNum decreases over time?
    paths(arma::span(0,CurrentPartNum-1), iter) = newS;
    iter++;
  }
  // arma::uword best = weights.index_max();
  // arma::urowvec bestPath = paths.row(best);
  return List::create(Named("paths") = paths,
                      Named("weights") = weights,
                      Named("LastStep") = iter++);
}


List pathStuffold(List pmats, arma::uvec path, arma::mat y){
    // Note: this function's smoother isn't quite right, see below
    // What if I want different initial state (instead of 0)?
    
    arma::mat a0 = pmats["a0"];
    arma::mat P0 = pmats["P0"];
    arma::cube dt = pmats["dt"];
    arma::cube ct = pmats["ct"];
    arma::cube Tt = pmats["Tt"];
    arma::cube Zt = pmats["Zt"];
    arma::cube Rt = pmats["Rt"];
    arma::cube Qt = pmats["Qt"];
    arma::cube GGt = pmats["GGt"];
    
    arma::uword n = y.n_cols;
    arma::uword m = a0.n_rows;
    arma::uword mm = m*m;
    arma::uword d = y.n_rows;
    arma::uword dm = d*m;
    arma::uword dd = d*d;
    
    arma::uword dtvar = dt.n_slices > 1;
    arma::uword ctvar = ct.n_slices > 1;
    arma::uword Ttvar = Tt.n_slices > 1;
    arma::uword Ztvar = Zt.n_slices > 1;
    arma::uword Rtvar = Rt.n_slices > 1;
    arma::uword Qtvar = Qt.n_slices > 1;
    arma::uword GGtvar = GGt.n_slices > 1;
    arma::mat HHt(m,m);
    arma::mat R(mm,1);
    arma::mat Q(mm,1);
    
    
    // output storage
    arma::colvec llik = arma::zeros(n);
    arma::mat at(m,n,arma::fill::zeros); // predicted mean E[a_t | y_1:t-1]
    arma::cube Pt(m,m,n,arma::fill::zeros); // predicted variance
    arma::mat preds(d,n,arma::fill::zeros); // predictions  
    arma::mat ahat(m,n,arma::fill::zeros); // smoothed mean E[a_t | y_1:n]
    arma::cube Phat(m,m,n,arma::fill::zeros); // smoothed variance
    arma::mat ests(d,n,arma::fill::zeros); //smoothed estimates
    double liktmp = 0.0;
    
    // initialization
    arma::uword s = path(0); // See comment above, replace 0 with appropriate initial s.
    arma::mat a00 =  a0.col(s);
    arma::mat P00 = reshape(P0.col(s), m, m);
    // arma::mat Z0 = Zt.subcube(0,s,0,arma::size(dm,1,1));
    // Z0.reshape(d,m);
    // arma::colvec c0 = ct.subcube(0,s,0,arma::size(d,1,1)); // does this work??
    // preds.col(0) = c0 + Z0 * at.col(0);
        
    
    // Kalman filtering
    for(arma::uword iter=0; iter<n; iter++){ //issue only happens when variables are declared both in and before the loop
        s = path(iter);
        if(iter==0 || Rtvar || Qtvar){ 
            R = Rt.subcube(0,s,iter*Rtvar,arma::size(mm,1,1));
            Q = Qt.subcube(0,s,iter*Qtvar,arma::size(mm,1,1));
            R.reshape(m,m);
            Q.reshape(m,m);
            HHt = R * Q * R.t();
        }
        
        KFOUT step = kf1step(a00, P00, 
                            dt.subcube(0,s,iter*dtvar,arma::size(m,1,1)),
                            ct.subcube(0,s,iter*ctvar,arma::size(d,1,1)), 
                            Tt.subcube(0,s,iter*Ttvar,arma::size(mm,1,1)),
                            Zt.subcube(0,s,iter*Ztvar,arma::size(dm,1,1)), 
                            HHt, GGt.subcube(0,s,iter*GGtvar,arma::size(dd,1,1)), 
                            y.col(iter));
        at.col(iter) = step.at;
        Pt.slice(iter) = arma::reshape(step.Pt, m, m);
        preds.col(iter) = step.pred;
        a00 = step.att;
        P00 = arma::reshape(step.Ptt, m, m);
        liktmp = step.lik;
        llik(iter) += log(liktmp);
    }
    
    // Kalman smoothing. Recalculation is inefficient, but likely doesn't matter.
    arma::mat r1(m,1,arma::fill::zeros);
    arma::mat N1(m,m,arma::fill::zeros);
    arma::uword iter = n;
    while(iter > 0){
        iter--;
        s = path(iter);
        arma::mat PP = Pt.slice(iter); // for ease
        arma::mat cc = ct.subcube(0,s,iter*ctvar,arma::size(d,1,1));
        arma::mat ZZ = Zt.subcube(0,s,iter*Ztvar,arma::size(dm,1,1));
        KFOUT step = ks1step(r1, N1, at.col(iter), PP, 
                             cc, Tt.subcube(0,s,iter*Ttvar,arma::size(mm,1,1)),
                             ZZ, GGt.subcube(0,s,iter*GGtvar,arma::size(dd,1,1)), 
                             y.col(iter));
        r1 = step.at;
        N1 = step.Pt;
        ahat.col(iter) = at.col(iter) + PP * r1;
        Phat.slice(iter) = PP - PP * N1 * PP;
        ests.col(iter) = cc + arma::reshape(ZZ, d, m) * ahat.col(iter);
    }
    return List::create(Named("at") = at, Named("Pt") = Pt, Named("preds") = preds,
                              Named("ahat") = ahat, Named("Phat") = Phat,
                              Named("ests") = ests, Named("llik") = llik);
}



//' Estimate continuous states given parameters and discrete hidden states
//' 
//' @param pmats e.g., as output from yupengMats
//' @param path path of discrete hidden states
//' @param y observations
//' 
//' @return List with components from Kalman filter and Smoother
//' 
//' @export 
// [[Rcpp::export]]
List pathStuff(List pmats, arma::uvec path, arma::mat y){
  // What if I want different initial state (instead of 0)?
  
  arma::mat a0 = pmats["a0"];
  arma::mat P0 = pmats["P0"];
  arma::cube dt = pmats["dt"];
  arma::cube ct = pmats["ct"];
  arma::cube Tt = pmats["Tt"];
  arma::cube Zt = pmats["Zt"];
  arma::cube Rt = pmats["Rt"];
  arma::cube Qt = pmats["Qt"];
  arma::cube GGt = pmats["GGt"];
  
  arma::uword n = y.n_cols;
  arma::uword m = a0.n_rows;
  arma::uword mm = m*m;
  arma::uword d = y.n_rows;
  arma::uword dm = d*m;
  arma::uword dd = d*d;
  
  arma::uword dtvar = dt.n_slices > 1;
  arma::uword ctvar = ct.n_slices > 1;
  arma::uword Ttvar = Tt.n_slices > 1;
  arma::uword Ztvar = Zt.n_slices > 1;
  arma::uword Rtvar = Rt.n_slices > 1;
  arma::uword Qtvar = Qt.n_slices > 1;
  arma::uword GGtvar = GGt.n_slices > 1;
  arma::mat HHt(m,m);
  arma::mat R(mm,1);
  arma::mat Q(mm,1);
  
  
  // output storage
  arma::colvec llik = arma::zeros(n);
  arma::mat at(m,n,arma::fill::zeros); // predicted mean E[a_t | y_1:t-1]
  arma::mat att(m,n,arma::fill::zeros); // updated mean E[a_t | y_1:t]
  arma::cube Pt(m,m,n,arma::fill::zeros); // predicted variance
  arma::cube Ptt(m,m,n,arma::fill::zeros); // predicted variance
  arma::mat preds(d,n,arma::fill::zeros); // predictions  
  arma::mat ahat(m,n,arma::fill::zeros); // smoothed mean E[a_t | y_1:n]
  arma::mat ests(d,n,arma::fill::zeros); //smoothed estimates
  double liktmp = 0.0;
  
  // initialization
  arma::uword s = path(0); // See comment above, replace 0 with appropriate initial s.
  arma::mat a00 =  a0.col(s);
  arma::mat P00 = reshape(P0.col(s), m, m);
  
  
  // Kalman filtering
  arma::uword iter=0;
  while(iter < n){ //issue only happens when variables are declared both in and before the loop
    s = path(iter);
    if(iter==0 || Rtvar || Qtvar){ 
      R = Rt.subcube(0,s,iter*Rtvar,arma::size(mm,1,1));
      Q = Qt.subcube(0,s,iter*Qtvar,arma::size(mm,1,1));
      R.reshape(m,m);
      Q.reshape(m,m);
      HHt = R * Q * R.t();
    }
    
    KFOUT step = kf1step(a00, P00, 
                         dt.subcube(0,s,iter*dtvar,arma::size(m,1,1)),
                         ct.subcube(0,s,iter*ctvar,arma::size(d,1,1)), 
                         Tt.subcube(0,s,iter*Ttvar,arma::size(mm,1,1)),
                         Zt.subcube(0,s,iter*Ztvar,arma::size(dm,1,1)), 
                         HHt, GGt.subcube(0,s,iter*GGtvar,arma::size(dd,1,1)), 
                         y.col(iter));
    at.col(iter) = step.at;
    att.col(iter) = step.att;
    Pt.slice(iter) = arma::reshape(step.Pt, m, m);
    Ptt.slice(iter) = arma::reshape(step.Ptt, m, m);
    preds.col(iter) = step.pred;
    a00 = step.att;
    P00 = arma::reshape(step.Ptt, m, m);
    liktmp = step.lik;
    llik(iter) += log(liktmp);
    iter++;
  }
  
  // Kalman smoothing. Based on Durbin and Koopman (4.70)
  iter--;
  ahat.col(iter) = att.col(iter);
  while(iter > 0){
    arma::mat cc = ct.subcube(0,s,iter*ctvar,arma::size(d,1,1));
    arma::mat zz = Zt.subcube(0,s,iter*Ztvar,arma::size(dm,1,1));
    zz.reshape(d,m);
    ests.col(iter) = cc + zz * ahat.col(iter);
    P00 = arma::pinv(Pt.slice(iter));
    a00 = ahat.col(iter) - at.col(iter);
    
    iter--; // This is important, need T after the increment rather than before.
    arma::mat Tmat = Tt.subcube(0,s,iter*Ttvar,arma::size(mm,1,1));
    s = path(iter);
    Tmat.reshape(m,m);
    ahat.col(iter) = att.col(iter) + Ptt.slice(iter) * Tmat * P00 * a00;
  }
  // iter = 0, need the estimates
  arma::mat cc = ct.subcube(0,s,iter*ctvar,arma::size(d,1,1));
  arma::mat zz = Zt.subcube(0,s,iter*Ztvar,arma::size(dm,1,1));
  zz.reshape(d,m);
  ests.col(iter) = cc + zz * ahat.col(iter);
    
  return List::create(Named("at") = at, Named("Pt") = Pt, Named("preds") = preds,
                            Named("ahat") = ahat,
                            Named("att") = att, Named("Ptt") = Ptt,
                            Named("ests") = ests, Named("llik") = llik);
}
