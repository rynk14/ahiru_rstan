//
// リターンの誤差分布をt分布にする
// y_t = exp(ht/2)ep_t, ep_t ~ t(nu, 0, 1)
// h_t+1 = mu + phi(h_t-mu)+eta_t, eta_t~N(0, sigma^2)

data {
  int T;
  vector[T] y;
}

parameters {
  real mu;
  real<lower=0> sigma_eta;
  real<lower=-1, upper=1> phi;
  real<lower=4> nu;
  vector[T] h;
}

transformed parameters {
  real<lower=-1, upper=1> phi_trans;
  real<lower=0> sigma_eta_sqrt;
  phi_trans = (1+phi)/2;
  sigma_eta_sqrt = sqrt(sigma_eta);
}

model {
  // 事前分布
  mu ~ normal(-10, 1);
  phi_trans ~ beta(20, 1.5);
  sigma_eta_sqrt ~ gamma(1, 5);
  nu ~ exponential(0.33333);
  
  // モデル
  h[1] ~ normal(mu, sigma_eta/sqrt(1-square(phi)));
  for (t in 2:T) {
    h[t] ~ normal(mu+phi*(h[t-1]-mu), sigma_eta);
    
  }
  for (t in 1:T) {
    y[t] ~ student_t(nu, 0, exp(h[t]/2));
  }
}



