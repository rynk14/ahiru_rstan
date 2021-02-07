//
// 多変量SVモデル
// 行列のテクニックとかはこのスライドがすごくわかりやすい
// https://www.slideshare.net/KojiKosugi/stanrch9
// 分散の事前分布は設定しない

data {
  int dim;
  int T;
  matrix[dim, T] y;
}

parameters {
  vector<lower=0>[dim] sigma;
  vector<lower=0, upper=1>[dim] phi_trans;
  matrix[dim, T+1] h;
}

transformed parameters {
  vector<lower=-1, upper=1>[dim] phi; 
  cov_matrix[dim] Sigma_eta;
  cov_matrix[dim] Sigma_0;
  vector[2*dim] y_and_vola[T+1];
  
  // 事前分布のための変換
  phi = 2*phi_trans-1;
  Sigma_eta = diag_matrix(square(sigma));
  // MSV元祖p.4参照、Sigma0
  Sigma_0 = Sigma_eta * diag_matrix(1 ./ (1-square(phi)));
  
  y_and_vola[1][1:dim] = rep_vector(0, dim);
  y_and_vola[1][(dim+1):(2*dim)] = h[,1];
  for (t in 2:(T+1)) {
    y_and_vola[t][1:dim] = y[, t-1];
    y_and_vola[t][(dim+1):(2*dim)] = h[, t];
  }
}

model {
  // 事前パラメータ
  vector[2*dim] mu;
  matrix[2*dim, 2*dim] Sigma;
  
  // 事前分布
  phi_trans ~ beta(20, 1.5);

  y_and_vola[1][dim+1] ~ normal(0, sigma[1]);
  y_and_vola[1][dim+2] ~ normal(0, sigma[2]);
  
  // サンプリング
  // modelブロックで=で書けるのはtransformed parametersの変数だけ。
  for (t in 2:(T+1)){
    mu[1:dim] = rep_vector(0, dim);
    mu[(dim+1):(dim*2)] = phi .* y_and_vola[t-1][(dim+1):(dim*2)];
    Sigma[1:dim, 1:dim] = diag_matrix(exp(y_and_vola[t-1][(dim+1):(dim*2)]));
    Sigma[(dim+1):(2*dim), (dim+1):(2*dim)] = Sigma_eta;
    
    // リターンとボラが無相関としてるので対角は０行列
    Sigma[(dim+1):(2*dim), 1:dim] = rep_matrix(0, dim, dim);
    Sigma[1:dim, (dim+1):(2*dim)] = Sigma[(dim+1):(2*dim), 1:dim];
    
    // モデル
    y_and_vola[t] ~ multi_normal(mu, Sigma);
  }
  
}


