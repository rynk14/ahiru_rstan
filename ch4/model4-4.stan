// p.39 の単回帰モデル

data{
  int N;
  real X[N];
  real Y[N];
}

parameters{
  real a; // 実数パラメータa
  real b; // 実数パラメータb
  real<lower = 0> sigma; // 非負実数パラメータσ
}

model{
  for (n in 1:N){
    Y[n] ~ normal(a+b*X[n], sigma);
  }
}
