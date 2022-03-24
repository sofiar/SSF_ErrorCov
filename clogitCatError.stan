// clogit.stan categorical error predictor

data {
  int<lower=0> N;
  int<lower=0> n_grp;
  int qgr[N];
  vector[N] x;
  vector[N] x2;
  int<lower=1,upper=2> x_cat[N];
  int<lower=0,upper=1> y[N];
  matrix[2,2] lQ;


}
parameters {
  vector[2] beta;
  vector[n_grp]alpha;
  
}

// 
 transformed parameters {
vector[2] log_q_x[N];
   for (i in 1:N) {
  log_q_x[i, 1] = bernoulli_logit_lpmf(y[i]|alpha[qgr[i]]+beta[1]*x2[i])+lQ[x_cat[i],1];//log(Q[1,x[i]+1]);
  log_q_x[i, 2] = bernoulli_logit_lpmf(y[i]|alpha[qgr[i]]+beta[1]*x2[i]+beta[2])+lQ[x_cat[i],2];//log(Q[2,x[i]+1]);

 }
}

model {
 beta~normal(0,2);
 
 
for (n in 1:N)
{

  target += log_sum_exp(log_q_x[n]);

}
  
}
