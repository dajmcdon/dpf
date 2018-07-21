functions {
  real stanloglike(row_vector y, vector lt, real sig2eps, 
    real mu1, real mu2, real mu3, 
    real sig2obs, real sig2acc, real sig2dec, real sig2stress, 
    real p11, real p12, real p22, real p31);
  real beam_lpdf(row_vector y, vector lt, real sig2eps, 
    real mu1, real mu2, real mu3, 
    real sig2obs, real sig2acc, real sig2dec, real sig2stress, 
    real p11, real p12, real p22, real p31){
      real llike;
      llike = stanloglike(y, lt, sig2eps, mu1, mu2, mu3, 
        sig2obs, sig2acc, sig2dec, sig2stress, p11, p12, p22, p31);
    return llike;
  }
} 
data {
  int<lower=0> N;
  row_vector[N] y;
  vector[N] lt;
}
transformed data{
  real samp_mean;
  samp_mean = mean(y);
}
parameters {
  real<lower=0> sig2eps;
  real<lower=0> mu1;
  real<lower=0> mu2;
  real<lower=0> mu3;
  real<lower=0> sig2obs;
  real<lower=0> sig2acc;
  real<lower=0> sig2dec;
  real<lower=0> sig2stress;
  simplex[3] p1;
  real<lower=0,upper=1> p22;
  real<lower=0,upper=1> p31;
}
transformed parameters {
  real p11 = p1[1];
  real p12 = p1[2];
} 
model {
  vector[3] alp;
  alp[1] = 7.0;
  alp[2] = 2.0;
  alp[3] = 1.0;
  sig2eps ~ normal(400.0, 100.0);
  mu1 ~ normal(samp_mean, 20.0);
  mu2 ~ normal(-40.0, 20.0);
  mu3 ~ normal(-40.0, 20.0);
  sig2obs ~ normal(0.0001, 0.1);
  sig2acc ~ normal(400, 100);
  sig2dec ~ normal(400, 100);
  sig2stress ~ normal(900, 100);
  p1 ~ dirichlet(alp);
  p22 ~ beta(9,6);
  p31 ~ beta(10,10);
  target += beam_lpdf(y | lt, sig2eps,
                mu1, mu2, mu3, 
                sig2obs, sig2acc, sig2dec, sig2stress,
                p11, p12, p22, p31);
}
