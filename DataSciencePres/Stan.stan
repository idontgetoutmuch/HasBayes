data {
  int<lower=0> N;
  real x[N];
}

parameters {
  real mu;
  real<lower=0,upper=1000> sigma;
}
model{
  x     ~ normal(mu, sigma);
  mu    ~ normal(0, 1000);
}
