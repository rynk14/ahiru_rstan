//
// GARCH(1,1)をやってみる

// 観測できるデータをdataブロック
data {
  int<lower=0> T; # 時点の数
  real r[T]; # 収益率
  real<lower=0> sigma1; # ボラの初期値
}


parameters {
  real mu; # 収益率の平均
  real<lower=0> omega; # GARCHの切片項
  real<lower=0,upper=1> alpha; # GARCHの係数
  real<lower=0,upper=(1-alpha)> beta; # GARCHの1時点前のボラの係数
}

// 上で宣言したパラメータを使って違うパラメータσを作る
transformed parameters {
  real<lower=0> sigma[T]; # ボラ
  sigma[1] = sigma1; # ボラの初期値を代入
  # GARCHのボラの式
  for (t in 2:T)
    sigma[t] = sqrt(omega
                     + alpha * pow(r[t-1] - mu, 2)
                     + beta * pow(sigma[t-1], 2));
}

model {
  # リターンのモデル
  r ~ normal(mu,sigma);
}
// モデルのところは観測できるやつだけを書いとく？
// この場合ボラはモデルはあるが、観測できないからここには書かない
// modelブロックにボラのモデルどう書くねん (分布がわからん) となって迷った
// ボラはパラメータのところに書く
// まぁ確かにボラはリターンの標準偏差やしパラメータと言われればそう

// じゃあなんで12-2のモデルでトレンドの項μをモデルのとこに書いてるのてなる。
// μは観測できない、yは実際の来場者数だから観測可能
// たしかにmodel12-2において、μは平均パラメータではあるが、モデルブロックにも書いてる
// ということはmodelのとこは、~分布() 　の形で書ける。
// 要は分布がわかってる決まってるものだけをmodelブロックに書いてるのかな。。

