data {
  int<lower = 0> n[2]; // number of tested persons
  int<lower = 0> y[2]; // number of positive persons
  int<lower = 0> n_Sp[2]; // number of true negative persons tested 
  int<lower = 0> y_Sp[2]; // number of true negative persons tested negative
  int<lower = 0> n_Se[2]; // number of true positive persons tested 
  int<lower = 0> y_Se[2]; // number of true positive persons tested positive
  vector[2] mu;
  vector<lower = 0>[2] prev_beta;
  cov_matrix[2] Sigma;
}
parameters {
  real<lower = 0, upper = 1> prev[2];
  vector[2] logit_Se_Sp_1;
  vector[2] logit_Se_Sp_2;
}
transformed parameters {
  real logit_Se_1 = logit_Se_Sp_1[1];
  real logit_Sp_1 = logit_Se_Sp_1[2];
  
  real logit_Se_2 = logit_Se_Sp_2[1];
  real logit_Sp_2 = logit_Se_Sp_2[2];
  
  real<lower=0,upper=1> Se[2];
  real<lower=0,upper=1> Sp[2];
  real<lower = 0, upper = 1> p[2];
  real<lower = -1, upper = 1> diff;
  
  Se[1] = inv_logit(logit_Se_1);
  Sp[1] = inv_logit(logit_Sp_1);

  Se[2] = inv_logit(logit_Se_2);
  Sp[2] = inv_logit(logit_Sp_2);
  
  p[1] = prev[1] * Se[1] + (1 - prev[1]) * (1 - Sp[1]);
  p[2] = prev[2] * Se[2] + (1 - prev[2]) * (1 - Sp[2]);
  
  diff = prev[1] - prev[2];
  
}
model {
  y[1] ~ binomial(n[1], p[1]);
  y_Sp[1] ~ binomial(n_Sp[1], Sp[1]);
  y_Se[1] ~ binomial(n_Se[1], Se[1]);
  logit_Se_Sp_1 ~ multi_normal(mu, Sigma);
  prev[1] ~ beta(prev_beta[1],prev_beta[2]);
  
  y[2] ~ binomial(n[2], p[2]);
  y_Sp[2] ~ binomial(n_Sp[2], Sp[2]);
  y_Se[2] ~ binomial(n_Se[2], Se[2]);
  logit_Se_Sp_2 ~ multi_normal(mu, Sigma);
  prev[2] ~ beta(prev_beta[1],prev_beta[2]);
}
