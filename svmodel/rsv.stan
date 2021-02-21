//
// ========Realized Stochastic Volatility model====================

// 単変量SVモデルを実装する。
// 参考は大森先生のreview論文

data {
  int T;
  vector[T] y; // リターンの観測値
  vector[T] x; // 対数実現ボラ
}

parameters {
  vector[T+1] h; // 潜在対数ボラ
  real<lower=0, upper=1> phi_trans;
  real mu;
  real xi;
  real<lower=0, upper=1> rho_trans;
  real<lower=0> sigma_eta_sq;
  real<lower=0> sigma_u_sq;
}

transformed parameters {
  matrix[3, (T+1)] yxh;
  real<lower=-1, upper=1> phi;
  real<lower=-1, upper=1> rho;
  real<lower=0> sigma_eta;
  real<lower=0> sigma_u;
  cov_matrix[3] S;
  
  rho = 2*rho_trans-1;
  phi = 2*phi_trans-1;
  sigma_eta = sqrt(sigma_eta_sq);
  sigma_u = sqrt(sigma_u_sq);
  yxh[1, 1] = 0;
  yxh[2, 1] = 0.5;
  yxh[1, 2:(T+1)] = y';
  yxh[2, 2:(T+1)] = x';
  yxh[3,] = h';
  
  S[1,1] = 1;
  S[1,2] = 0;
  S[1,3] = rho*sigma_eta;
  S[2,1] = 0;
  S[2,2] = sigma_u_sq;
  S[2,3] = 0;
  S[3,1] = rho*sigma_eta;
  S[3,2] = 0;
  S[3,3] = sigma_eta_sq;

}

model {
  vector[3] ave;
  matrix[3, 3] Sigma;
  vector[3] V;
  
  mu ~ normal(0, 10);
  phi_trans ~ beta(20, 1.5);
  sigma_eta_sq ~ inv_gamma(2.5, 0.025);
  rho_trans ~ beta(1, 2);
  xi ~ normal(0, 1);
  sigma_u_sq ~ inv_gamma(2.5, 0.1);
  
  h[1] ~ normal(mu, sigma_eta/sqrt(1-square(phi)));
  for (t in 2:(T+1)) {
    ave[1] = 0;
    ave[2] = xi+yxh[3,t-1];
    ave[3] = mu+phi*(h[t-1]-mu);
    V[1] = exp(h[t-1]/2);
    V[2:3] = rep_vector(1, 2);
    Sigma = diag_matrix(V) * S * diag_matrix(V);
    //print(Sigma);
    
    yxh[,t] ~ multi_normal(ave, Sigma);
  }
}

