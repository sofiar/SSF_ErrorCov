
data{
  int N; 
  vector[N] y;
  }

parameters{
  real<lower=0> beta;
}

model
{
 beta~normal(0,1);
  y ~ rayleigh(beta);
}

