data {
  int<lower=0> N1;
  int<lower=0> N2;
  real time_int[N1 + N2];
  
  int<lower=0> max_k;
  int<lower=0, upper=max_k> k[N1];
  matrix<lower=0>[N1, max_k] m;
  int<lower=0, upper=1> y[N1, max_k];
  
  real<lower=0> ig_alpha;
  real<lower=0> ig_beta;
}

transformed data {
  real delta = 1e-9;
  int<lower=1> N = N1 + N2;
}

parameters {
  real<lower=0> elleq;
  real<lower=0> sigma;
  vector[N] eta;
  real mu;
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
  for (i in 1:N1){
    for (j in 1:k[i]){
      y[i, j] ~ bernoulli(1 - (1 - Phi(z[i]))^m[i, j]); 
    }
  }
}

