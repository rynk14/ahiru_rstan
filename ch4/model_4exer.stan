// p.53 練習問題

data {
  int N1;
  int N2;
  real Y1[N1];
  real Y2[N2];
}

parameters {
  real mu1;
  real mu2;
  real<lower=0> sigma;
}

model {
  for (i in 1:N1)
    Y1[i] ~ normal(mu1, sigma);
  for (i in 1:N2)
    Y2[i] ~ normal(mu2, sigma);
}
