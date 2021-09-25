data {
  int<lower=0> N;
  real time_int[N];
  
  int<lower=0> k[N];
  int<lower=1> m_all;
  int<lower=0> m[N];
  
  int<lower=0> y[N, 2];
  
  real<lower=0> ig_alpha;
  real<lower=0> ig_beta;
}

transformed data {
  real delta = 1e-9;
}

parameters {
  real<lower=0> elleq;
  real<lower=0> sigma;
  real mu;
  vector[N] eta;
}

transformed parameters {
  vector[N] z;
  {
    matrix[N, N] L_K;
    matrix[N, N] K = cov_exp_quad(time_int, sigma, elleq);

    for (i in 1:N)
      K[i, i] = K[i, i] + delta;

    L_K = cholesky_decompose(K);
    z = L_K * eta + mu;
  }
}

model {
  elleq ~ inv_gamma(ig_alpha, ig_beta);
  sigma ~ normal(0, 0.5);
  eta ~ std_normal();
  mu ~ normal(-2, 0.5);
  for (i in 1:N){
    y[i,1] ~ binomial(k[i], 1 - (1 - Phi(z[i]))^m_all);
    y[i,2] ~ bernoulli(1 - (1 - Phi(z[i]))^m[i]);
  }
}
