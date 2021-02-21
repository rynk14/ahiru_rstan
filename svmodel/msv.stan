// なんか動いた
// Sigma_oriで共分散行列にしたのがよかった


data {
  int dim;
  int T;
  matrix[dim, T] y;
}

parameters {
  vector<lower=0, upper=5>[dim] sigma_eps;
  vector<lower=0, upper=5>[dim] sigma_eta;
  vector<lower=0, upper=1>[dim] phi_trans;
  matrix[dim, T+1] h;
  matrix<lower=-1, upper=1>[dim, dim] rho_etaeps; // i,j成分はi番目の対数ボラとj番目のリターンの相関係数
  // rho_etaepsは相関行列ではないことに注意
  corr_matrix[dim] rho_eta; // i,j成分はi番目の対数ボラとj番目の対数ボラの相関係数
  corr_matrix[dim] rho_eps; // i,j成分はi番目のリターンとj番目のリターンの相関係数
  cov_matrix[2*dim] Sigma_ori; // ベースになる共分散をここで宣言しても結果は同じ
}

transformed parameters {
  //cov_matrix[dim] Sigma_eps;
  //cov_matrix[dim] Sigma_eta; 
  //matrix[dim, dim] Sigma_etaeps; // こいつ自体は共分散行列ではない
  //cov_matrix[dim*2] Sigma_ori; // 共分散行列の時間に関係無い部分
  cov_matrix[dim] Sigma_0;
  vector[dim] mu_0;
  vector<lower=-1, upper=1>[dim] phi;
  vector[2*dim] y_and_vola[T+1];
  matrix[(2*dim), (2*dim)] Sigma_star;
  matrix[dim, dim] I_d ;
  matrix[dim, dim] I1;
  
  // 事前分布のための変換
  phi = 2*phi_trans-1;
  //Sigma_eps = (sigma_eps*(sigma_eps')) .* rho_eps;
  //Sigma_eta = (sigma_eta*(sigma_eta')) .* rho_eta;
  //Sigma_etaeps = (sigma_eps*(sigma_eta')) .* rho_etaeps;
  //Sigma_ori[1:dim, 1:dim] = Sigma_eps;
  //Sigma_ori[1:dim, (1+dim):(2*dim)] = Sigma_etaeps;
  //Sigma_ori[(1+dim):(2*dim), 1:dim] = Sigma_ori[1:dim, (1+dim):(2*dim)]';
  //Sigma_ori[(1+dim):(2*dim), (1+dim):(2*dim)] = Sigma_eta;
  
  //Sigma_ori[1:dim, 1:dim] = (sigma_eps*(sigma_eps')) .* rho_eps;
  //Sigma_ori[1:dim, (1+dim):(2*dim)] = (sigma_eps*(sigma_eta')) .* rho_etaeps;
  //Sigma_ori[(1+dim):(2*dim), 1:dim] = Sigma_ori[1:dim, (1+dim):(2*dim)]';
  //Sigma_ori[(1+dim):(2*dim), (1+dim):(2*dim)] = (sigma_eta*(sigma_eta')) .* rho_eta;
  //print(eigenvalues_sym(Sigma_ori));

  // MSV元祖p.4参照、Sigma0
  Sigma_0 = ((sigma_eta*(sigma_eta')) .* rho_eta) ./ (phi*(phi'));
  mu_0 = rep_vector(0, dim);
  
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
  vector[dim*2] mu;
  matrix[2*dim, 2*dim] Sigma;
  vector[2*dim] V;
  V[(1+dim):2*dim] = rep_vector(1, dim);
  // 事前分布
  phi_trans ~ beta(20, 1.5);
  
  // ボラの初期分布
  y_and_vola[1][(dim+1):(2*dim)] ~ multi_normal(mu_0, Sigma_0);
  
  // サンプリング
  // modelブロックで=で書けるのはmodelで定義したやつだけ。
  for (t in 2:(T+1)){
    mu[1:dim] = rep_vector(0, dim);
    mu[(dim+1):(dim*2)] = phi .* h[,t-1];
    V[1:dim] = exp(h[,t-1]/2);
    Sigma = diag_matrix(V) * Sigma_ori *diag_matrix(V);
    
    // モデル
    Sigma ~ inv_wishart(10, Sigma_star); ;// 共分散行列の事前分布
    y_and_vola[t] ~ multi_normal(mu, Sigma);
  }
}

