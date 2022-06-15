data {
  int<lower=0> N1;
  int<lower=0> N2;
  real time_int[N1 + N2];
  
  int<lower=0> k;
  int<lower=0> y[N1];
  
  int<lower=0> num_group;
  int<lower=0> group[N1];

  real<lower=0> ig_alpha;
  real<lower=0> ig_beta;
}

transformed data {
  real delta = 1e-9;
  int<lower=1> N = N1 + N2;
  int<lower=N1> N_start = N1 + 1;
}

parameters {
  real<lower=0> elleq;
  real<lower=0> sigma;
  real mu;
  vector[N] eta;
  vector[num_group] mu_group;
}

transformed parameters {
  vector[N] z;
  vector[N] mu_vec;
  
  for (i in N_start:N){
    mu_vec[i] = mu;
  }
  
  for (i in 1:N1){
    mu_vec[i] = mu_group[group[i]];
  }

  {
    matrix[N, N] L_K;
    matrix[N, N] K = cov_exp_quad(time_int, sigma, elleq);

    for (i in 1:N)
      K[i, i] = K[i, i] + delta;

    L_K = cholesky_decompose(K);
    z = L_K * eta + mu_vec;
  }
}

model {
  elleq ~ inv_gamma(ig_alpha, ig_beta);
  sigma ~ normal(0, 0.5);
  eta ~ std_normal();
  mu ~ normal(-2, 0.5);
  y ~ binomial(k, Phi(z[1:N1]));
  mu_group ~ normal(mu, .2);
}
