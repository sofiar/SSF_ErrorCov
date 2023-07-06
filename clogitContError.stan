
data {
  int<lower=0> N; // Number of observations
  int<lower=1> n_grp; // Number of groups
  int qgr[N];
  int<lower=1> n_coef; // Number of coefficients (log odds ratios) to estimate
  int<lower=0, upper=1> y[N]; // array of 0/1 outcomes
  matrix[N, n_coef] x; // Matrix of regressors
  real<lower=0> tau;  // measurement noise
  real<lower=0> sigma_x;  // sd x
  
}

parameters {
  vector[n_coef] b;
  vector[n_grp]alpha;
  vector[N] xt;          // unknown true values
 }

model {
  //vector[N] xb; // observation level linear predictor
  vector[N] xm;
  vector[N] x2;
  
  xm = col(x,1);
  x2 = col(x,2);
  
  
  
  alpha~normal(-5,10);
  // diffuse normal prior for log odds ratios
  b ~ normal(0, 3);
  //error measurement
  xt ~ normal(0, sigma_x);  // prior
  xm ~ normal(xt, tau);    // measurement model
  
  //xb = xt*b[1]+x2*b[2];
  
for (i in 1:N)
{

  //y[i] ~ bernoulli_logit(alpha[qgr[i]] + b[1]*xt[i] +b[2]*x2[i]);
  target += bernoulli_logit_lpmf(y[i]|alpha[qgr[i]]+xt[i]*b[1]+x2[i]*b[2]);

}








  // log likelihood is a sum over each group
  //pos = 1;
  //for (ii in 1:n_grp) {
  //  int y_g[n_group[ii]];
  //  vector[n_group[ii]] xb_g;
  //  y_g = segment(y, pos, n_group[ii]);
  //  xb_g = segment(xb, pos, n_group[ii]);
  //  ll = dot_product(to_vector(y_g), xb_g+alpha[ii]) - log(cl_denom(n_group[ii], n_case[ii], xb_g+alpha[ii]));
  //  target += ll;
  //  pos = pos + n_group[ii];
  //}
}

