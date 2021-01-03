// p.50 model4-4.stanの写す

data {
  int N;
  real X[N];
  real Y[N];
  int N_new;
  real X_new[N_new];
}

parameters {
  real a;
  real b;
  real<lower=0> sigma;
}

// parameterとdataを使って新たな変数を作る
transformed parameters {
  real y_base[N];
  for (n in 1:N){
    y_base[n] = a + b*X[n];
  }
}

model {
  for (n in 1:N){
    Y[n] ~ normal(y_base[n], sigma);
  }
}

// data、parameters, transformed parametersブロックで宣言されたやつらを使って
// 新しい変数を定義する。このブロックは事後確率とは完全に切り離されていて計算が速い
// このブロック内で従うを意味する ~ がつかえないから、ある分布に従う乱数を発生
// させる場合は、~分布名() の代わりに、=分布名_rng()をつかう
// ここでは新しいデータのやつを宣言
generated quantities {
  real y_base_new[N_new];
  real y_new[N_new];
  for (n in 1:N_new){
    y_base_new[n] = a + b*X_new[n];
    y_new[n] = normal_rng(y_base_new[n], sigma);
  }
  
}

