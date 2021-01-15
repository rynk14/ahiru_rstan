//　今度はSVモデルもやってみる
// ただSVモデルは これ！ といった決まりがなくて色んなモデルがある。
// ここでは大森先生のhttps://www.terrapub.co.jp/journals/jjssj/pdf/4202/42020273.pdf
// に載ってるやつを実装する
// 簡単のために二つの誤差項同士の相関ρは0にした


// こんな風にベクトルを駆使すると高速になるらしい

data {
  int T; // 観測値の数, 時点の数でもある
  vector[T] ret; // 対数リターンのベクトル
}

parameters {
  real mu; // 平均
  real<lower=0> sigma; // 誤差の標準偏差
  vector[T] h; // 対数ボラ
  real<lower=-1, upper=1> phi; // 持続性パラメータ
}

model {
  h[1] ~ normal(mu, sigma/sqrt(1-phi*phi));
  h[2:T] ~ normal(mu+phi*(h[1:(T-1)]-mu), sigma);
  ret ~ normal(0, exp(h/2));
}

