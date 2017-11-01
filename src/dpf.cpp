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
    arma::mat a1;
    arma::mat P1;
    double lik;
    arma::mat pred;
};

KFOUT kf1step(arma::mat a0, arma::mat P0, arma::mat dt,
             arma::mat ct, arma::mat Tt,
             arma::mat Zt, arma::mat HHt, arma::mat GGt, arma::mat yt) {
  // Reshape matrices
  int d = a0.size();//number of elements.
  int m = yt.size();
  P0.reshape(d,d);
  Tt.reshape(d,d);
  Zt.reshape(m,d);
  HHt.reshape(d,d);
  GGt.reshape(m,m);

  // KF
  arma::mat a1 = a0;
  //Rcout << "line 138:" << a1 << std::endl;
  arma::mat pred = ct + Zt * a0;
  arma::mat vt = yt - pred;
  arma::mat Ft = GGt + Zt * P0 * Zt.t();
  arma::mat Ftinv = arma::inv(Ft);
  arma::mat Kt = P0 * Zt.t() * Ftinv;
  a1 += Kt * vt;
  //Rcout << "line 145:" << a1 << std::endl;
  arma::mat P1 = P0 - P0 * Zt.t() * Kt.t();
  a1 = dt + Tt * a1;
  //Rcout << "line 148:" << a1 << std::endl;
  P1 = HHt + Tt * P1 * Tt.t();

  // Calculate likelihood
  double ftdet = arma::det(Ftinv);
  ftdet = 1 / std::abs(ftdet);
  double mahalanobis = arma::as_scalar(vt.t() * Ftinv * vt);
  mahalanobis += yt.size()* log(2*M_PI) + log(ftdet);
  double lik = exp(-1.0*mahalanobis/2);
  P1.reshape(d*d, 1);
  KFOUT output = {a1, P1, lik, pred};
  return output;
  /*return List::create(Named("a1") = a1,
                      Named("P1") = P1,
                      Named("lik") = lik,
                      Named("pred") = pred);*/
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
//[[Rcpp::export]]
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
      arma::mat atmp = kfout.a1;
      arma::mat Ptmp = kfout.P1;
      a11.col(k) = atmp;
      P11.col(k) = Ptmp;
      lik(part,k) = kfout.lik;
    }
    a1.slice(part) = a11;
    P1.slice(part) = P11;
  }
  // arma::mat blah = transProbs.rows(currentStates);
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
        oldstates(i) = part;                                    //Shouldn't this be 'currentStates(part)'
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
                      Named("newW") = newW);
}



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
  
  arma::mat aa0 = a0.col(0);
  arma::mat PP0 = reshape(P0.col(0), m, m);
  
  double llik = 0;
  
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
    Rcout << step.a1 << std::endl;
    /*Rcout << "iter = " << iter << std::endl;
    Rcout << "s = " << s << std::endl;
    Rcout << "aa0 = " << aa0 << std::endl;
    Rcout << "PP0 = " << PP0 << std::endl;
    Rcout << "ct = " << ct.subcube(0,s,iter*ctvar,arma::size(d,1,1)) << std::endl;
    Rcout << "Tt = " << Tt.subcube(0,s,iter*Ttvar,arma::size(mm,1,1)) << std::endl;
    Rcout << "Zt = " << Zt.subcube(0,s,iter*Ztvar,arma::size(dm,1,1)) << std::endl;
    Rcout << "HHt = " << HHt << std::endl;
    Rcout << "GGt = " << GGt.subcube(0,s,iter*GGtvar,arma::size(dd,1,1)) << std::endl;
    Rcout << "yt = " << y.col(iter) << std::endl;*/
    
    arma::mat aa0 = step.a1;
    arma::mat PP0 = step.P1;
    double liktmp = step.lik;
    llik -= log(liktmp);
  }
  return llik;
}
    

    
//mus(0) = mu1, mus(1) = B section tempo, mus(2) = taot, mus(3) = phit
    
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
  dt.tube(0,4) = lt*mus(1);
  dt.tube(1,1) += mus(1);
  dt.tube(1,4) += mus(1);
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
  Tt.tube(2,4) = lt;
  Tt.tube(2,6) = lt;
  Tt.tube(2,3) = lt;
  arma::cube Zt(m, nstates, 1, arma::fill::zeros);
  Zt.tube(0,0,0,7) += 1;
  Zt(1,2,0) += 1;
  arma::cube Rt(mm, nstates, n, arma::fill::zeros);
  Rt.tube(0,5) += 1;
  Rt.tube(3,1,3,2) += 1;
  Rt.tube(3,4) += 1;
  Rt.tube(2,1) = lt;
  Rt.tube(2,4) = lt;
  arma::cube Qt(mm, nstates, 1, arma::fill::zeros);
  Qt.tube(0,0,0,7) += sig2eta(0);
  Qt.tube(3,0,3,7) += sig2eta(1);
  Qt.tube(3,2) -= sig2eta(1) - sig2eta(2);
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
  particles.head(K) = filler;
  arma::umat paths(maxpart, n, arma::fill::zeros);
  paths(arma::span(0,K-1),0) = filler;
  arma::colvec weights(maxpart,arma::fill::zeros);
  for(arma::uword iter=0; iter<K; iter++) weights(iter) += w0(iter);
  // weights.resize(maxpart);
  arma::uword CurrentPartNum = arma::accu(weights != 0);
  arma::mat HHt(m*m, K);
  arma::uword iter = 0;
  while(iter < n){
    if(iter==0 || Rtvar || Qtvar){ 
      HHt = HHcreate(Rt.slice(iter * Rtvar), Qt.slice(iter * Qtvar), m, q);
    }
    List step = dpf(particles.head(CurrentPartNum), weights.head(CurrentPartNum), 
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
    paths.head_rows(CurrentPartNum) = tmp;
    paths(arma::span(0,CurrentPartNum-1), iter) = newS;
    iter++;
  }
  // arma::uword best = weights.index_max();
  // arma::urowvec bestPath = paths.row(best);
  return List::create(Named("paths") = paths,
                      Named("weights") = weights,
                      Named("LastStep") = iter++);
}

// [[Rcpp::export]]
List pathStuff(List pmats, arma::uvec path, arma::mat y){
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
    
    arma::mat aa0 = a0.col(0);
    arma::mat PP0 = reshape(P0.col(0), m, m);
    
    arma::colvec llik = arma::zeros(n);
    
    arma::colvec means = arma::zeros(n);
    
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
        means(iter) = step.pred(0,0);
        arma::mat aa0 = step.a1;
        arma::mat PP0 = step.P1;
        double liktmp = step.lik;
        llik(iter) += log(liktmp);
    }
    return List::create(Named("preds") = means,
                        Named("llik") = llik);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/***R
# nstates = 2
#npart = 8
# N = 35
# n = 100
# y = matrix(1:3/10)
# y = matrix(rnorm(3*n),3)
# y = matrix(rnorm(n)+25,1)
# a0 = matrix(rnorm(2*npart),2)
# P0 = matrix(c(diag(.1,2)),4,npart)
# w0 = runif(npart)
# w0 = w0/sum(w0)
# ct = array(0,c(3,2,1))
# dt = array(0,c(2,2,1))
# Zt = array(c(.1,.1,.1,.2,.3,.1),c(6,2,1))
# Tt = array(c(diag(.9,2),diag(.8,2)),c(4,2,1))
# Rt = array(c(diag(1,2)),c(4,2,1))
# Qt = array(c(diag(1,2)),c(4,2,1))
# GGt = array(c(diag(.1,3)),c(9,2,1))
# currentStates = rbinom(npart,size = nstates-1,.5)
# w = runif(npart)
# w[2] = npart
#w = w/sum(w)
# currentStates[2] = 0
# transProbs = matrix(c(.1,.9,.4,.6),2,byrow=TRUE)
# out1 = dpf(currentStates, w, N, transProbs, a0, P0, dt, ct, Tt, Zt, HHt, GGt, y)
# ptm = proc.time()
# out = beamSearch(a0, P0, dt, ct, Tt, Zt, Rt, Qt, GGt, y, transProbs, N)
# print(proc.time() - ptm)
# lt = runif(n)
# temposwitch = double(n)
# temposwitch[floor(n/2):floor(3*n/4)] = 1
mus = c(60, 100, 2, 1)
# sig2eps = 1
# sig2eta = c(2, .1, 1)
# transProbs = c(.8,.1,.5,.4)
# ptm = proc.time()
# out = yupengMats(lt, temposwitch, sig2eps, mus, sig2eta, transProbs)
# test = beamSearch(a0,P0,w0, out$dt, out$ct, out$Tt, out$Zt, out$Rt, 
#                  out$Qt, out$GGt, y, out$transMat, 50)
# print(proc.time() - ptm)
toab <- function(x, a, b) x*(b-a) + a # maps [0,1] to [a,b]
logistic <- function(x) 1/(1+exp(-x)) # maps R to [0,1]


toOptimize <- function(pvec, lt, temposwitch, y, w0, Npart){
  sig2eps = exp(pvec[1])
  mus = pvec[2:5]
  sig2etas = exp(pvec[6:8])
  transprobs = logistic(pvec[9:12])
  transprobs[2] = toab(transprobs[2], 0, 1-transprobs[1])
  pmats = yupengMats(lt, temposwitch, sig2eps, mus, sig2etas, transprobs)
  S = beamSearch(pmats$a0, pmats$P0, w0, pmats$dt, pmats$ct, pmats$Tt, pmats$Zt,
                 pmats$Rt, pmats$Qt, pmats$GGt, y, pmats$transMat, Npart)
  if(S$LastStep < ncol(y)) return(Inf)
  best = S$paths[which.max(S$weights),]
  negllike = getloglike(pmats, best, y)
  return(negllike)
}
*/

