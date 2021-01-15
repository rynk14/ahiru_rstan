
// p.232のmodel12-2

data {
  int T;
  int T_pred;
  vector[T] Y;
}

parameters {
  vector[T] mu;
  real<lower=0> s_mu;
  real<lower=0> s_y;
}

model {
  mu[2:T] ~ normal(mu[1:(T-1)], s_mu);
  Y ~ normal(mu, s_y);
}

// 予測データはこのブロックに書く
generated quantities {
  vector[T+T_pred] mu_all; // 過去と予測分のmu
  vector[T_pred] y_pred; // 予測分の測定値
  mu_all[1:T] = mu; // 過去のデータ分はそのまま入れとく
  for (t in 1:T_pred) {
    mu_all[T+t] = normal_rng(mu_all[T+t-1], s_mu);
    y_pred[t] = normal_rng(mu_all[T+t], s_y);
  }
}