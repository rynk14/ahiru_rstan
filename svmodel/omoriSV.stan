// 大森先生の論文で用いられた単純な単変量ＳＶモデルを実装してみる。
// y_t = exp(h_t/2)ε_t
// h_t+1 = μ+φ(h_t - μ)+η_t
// h_1 ~ N(μ, σ_η^2/(1-φ^2))
// (ε_t η_t) ~ N_2{(0 0), ([1 ρσ_η], [ρσ_η σ_η^2])}

data {
  int T; // サンプルサイズ
  vector[T] y; // リターンの観測ベクトル
}

parameters {
  real mu_vola;
  real<lower=-1, upper=1> phi;
  real<lower=0> sigma_eta;
  real<lower=-1, upper=1> rho;
  vector[T] h;
}

transformed parameters {
  real<lower=-1, upper=1> phi_trans;
  real<lower=0> sigma_eta_square;
  real<lower=0, upper=1> rho_trans;
  vector[2] mu[T];
  cov_matrix[2] Sigma[T];
  matrix[T, 2] y_and_h;
  
  // 事前分布で使う変換
  phi_trans = (1+phi)/2;
  sigma_eta_square = square(sigma_eta);
  rho_trans = (1+rho)/2;
  
  mu[1,1] = 0;
  mu[1,2] = mu_vola;
  Sigma[1,1,1] = exp(h[1]);
  Sigma[1,1,2] = rho*exp(h[1]/2)*sigma_eta;
  Sigma[1,2,1] = Sigma[1,1,2];
  Sigma[1,2,2] = sigma_eta_square;
  for (t in 2:T) {
    mu[t,1] = 0;
    mu[t,2] = mu_vola+phi*(h[t-1]-mu_vola);
    Sigma[t,1,1] = exp(h[t-1]);
    Sigma[t,1,2] = rho*exp(h[t-1]/2)*sigma_eta;
    Sigma[t,2,1] = Sigma[t,1,2];
    Sigma[t,2,2] = sigma_eta_square;
  }
  
  y_and_h = append_col(y, h);
}

model {
  // 事前分布
  mu_vola ~ normal(0.4, 1);
  phi_trans ~ beta(20, 1.5);
  sigma_eta_square ~ inv_gamma(2.5, 0.025);
  rho_trans ~ beta(1, 2);
  
  // サンプリング
  y_and_h[1,] ~ multi_normal(mu[1,], Sigma[1,]);
  for (t in 2:T) {
    y_and_h[t,] ~ multi_normal(mu[t,], Sigma[t,]);
  }
}



