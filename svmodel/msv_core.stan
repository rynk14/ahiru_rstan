//
//
// 多変量SVモデル
// 行列のテクニックとかはこのスライドがすごくわかりやすい
// https://www.slideshare.net/KojiKosugi/stanrch9
// Sigmaの事前分布を設定した三田
// リターンとボラの相関を考慮する。
// 具体的には分散共分散行列の非対角を0以外にする
// これが第一目標の標準的なMSVモデルのつもり


data {
  int dim;
  int T;
  matrix[dim, T] y;
}

parameters {
  vector<lower=0>[dim] sigma;
  vector<lower=0, upper=1>[dim] phi_trans;
  matrix[dim, T+1] h;
  matrix<lower=-1, upper=1>[dim, dim] rho; // i,j成分はi番目のボラとj番目のリターンの相関係数
}

transformed parameters {
  vector<lower=-1, upper=1>[dim] phi; 
  cov_matrix[dim] Sigma_eta;
  matrix[dim, dim] Sigma_eta_sqrt;
  cov_matrix[dim] Sigma_0;
  vector[2*dim] y_and_vola[T+1];
  matrix[(2*dim), (2*dim)] Sigma_star;
  matrix[dim, dim] I_d ;
  matrix[dim, dim] I1;
  
  
  // 事前分布のための変換
  phi = 2*phi_trans-1;
  Sigma_eta = diag_matrix(square(sigma));
  Sigma_eta_sqrt = diag_matrix(sigma);
  // MSV元祖p.4参照、Sigma0
  Sigma_0 = Sigma_eta * diag_matrix(1 ./ (1-square(phi)));
  y_and_vola[1][1:dim] = rep_vector(0, dim);
  y_and_vola[1][(dim+1):(2*dim)] = h[,1];
  for (t in 2:(T+1)) {
    y_and_vola[t][1:dim] = y[, t-1];
    y_and_vola[t][(dim+1):(2*dim)] = h[, t];
  }
  
  // Sigmaの事前分布のためのやつ
  I_d = diag_matrix(rep_vector(1, dim));
  I1 = rep_matrix(1, dim, dim);
  Sigma_star[1:dim, 1:dim] = square(1.5) * (0.5*I_d+0.5*I1);
  Sigma_star[(1+dim):(2*dim), 1:dim] = 1.5*0.2*(-0.1)*I_d;
  Sigma_star[1:dim, (1+dim):(2*dim)] = (Sigma_star[(1+dim):(2*dim), 1:dim])';
  Sigma_star[(1+dim):(2*dim), (1+dim):(2*dim)] = square(0.2)*(0.5*I_d+0.5*I1);
  Sigma_star = 10*Sigma_star;
}

model {
  // 事前パラメータ
  vector[2*dim] mu;
  matrix[2*dim, 2*dim] Sigma;
  
  // 事前分布
  phi_trans ~ beta(20, 1.5);
  
  // ボラの初期分布
  y_and_vola[1][(dim+1):(2*dim)] ~ normal(0, sigma);
  
  // サンプリング
  // modelブロックで=で書けるのはtransformed parametersの変数だけ。
  for (t in 2:(T+1)){
    mu[1:dim] = rep_vector(0, dim);
    mu[(dim+1):(dim*2)] = phi .* y_and_vola[t-1][(dim+1):(dim*2)];
    Sigma[1:dim, 1:dim] = diag_matrix(exp(y_and_vola[t-1][(dim+1):(dim*2)]));
    Sigma[(dim+1):(2*dim), (dim+1):(2*dim)] = Sigma_eta;
    
    // リターンとボラが相関ありとしてるので非対角に適当な行列を設定する
    Sigma[1:dim, (dim+1):(2*dim)] = diag_matrix(exp(y_and_vola[t-1][(dim+1):(2*dim)]/2))*rho*Sigma_eta_sqrt;
    Sigma[(dim+1):(2*dim), 1:dim] = (Sigma[1:dim, (dim+1):(2*dim)])';
    // モデル
    Sigma ~ inv_wishart(10, Sigma_star); ;// 共分散行列の事前分布
    y_and_vola[t] ~ multi_normal(mu, Sigma);
  }
  
}

