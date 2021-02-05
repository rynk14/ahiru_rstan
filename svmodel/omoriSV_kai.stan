// 大森先生の論文で用いられた単純な単変量ＳＶモデルを実装してみる。
// y_t = exp(h_t/2)ε_t
// h_t+1 = φh_t+η_t
// h_1 ~ N(0, σ_η^2/(1-φ^2))
// (ε_t η_t) ~ N_2{(0 0), ([1 ρσ_η], [ρσ_η σ_η^2])}

// 計算速度が遅いのでどうにか改善を試みる
// https://www.slideshare.net/simizu706/stan-62042940 参考
// どうやら相関を考慮する多変量正規分布の処理は遅い
// そこでウィシャート分布を使う。
// 正規分布に従う値の2乗和の分布がカイ二乗
// χ二乗分布の多変量版がウィシャート
// データは偏差積和行列 (共分散にサンプルサイズかけたやつ)
// てやろうと思ったけど、Wishartでやると、平均をサンプリングできない。
// だから、平均が無いモデルを考える。


data {
  int T; // サンプルサイズ
  vector[T] y; // リターンの観測ベクトル
}

parameters {
  real<lower=-1, upper=1> phi;
  real<lower=0> sigma_eta;
  real<lower=-1, upper=1> rho;
  vector[T] h;
}

transformed parameters {
  real<lower=-1, upper=1> phi_trans;
  real<lower=0> sigma_eta_square;
  real<lower=0, upper=1> rho_trans;
  cov_matrix[2] Sigma[T];
  matrix[2, 2] spd; // 偏差積和行列
  
  // 事前分布で使う変換
  phi_trans = (1+phi)/2;
  sigma_eta_square = square(sigma_eta);
  rho_trans = (1+rho)/2;
  
  Sigma[1,1,1] = exp(h[1]);
  Sigma[1,1,2] = rho*exp(h[1]/2)*sigma_eta;
  Sigma[1,2,1] = Sigma[1,1,2];
  Sigma[1,2,2] = sigma_eta_square;
  for (t in 2:T) {
    Sigma[t,1,1] = exp(h[t-1]);
    Sigma[t,1,2] = rho*exp(h[t-1]/2)*sigma_eta;
    Sigma[t,2,1] = Sigma[t,1,2];
    Sigma[t,2,2] = sigma_eta_square;
  } 
  
  // 偏差積和行列を作る
  spd[1,1] = (T-1)*variance(h);
  spd[1,2] = sum((h-mean(h)).*(y-mean(y))); //.*で要素ごとの積
  spd[2,1] = spd[1,2];
  spd[2,2] = (T-1)*variance(y);

}

model {
  // 事前分布
  phi_trans ~ beta(20, 1.5);
  sigma_eta_square ~ inv_gamma(2.5, 0.025);
  rho_trans ~ beta(1, 2);
  
  // サンプリング
  for (t in 1:T) {
    spd ~ wishart(T-1, Sigma[t]); // Wishartはサンプルサイズと分散
  }
}


