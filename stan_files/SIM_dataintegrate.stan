data {
  int<lower=0> N_indiv;
  int<lower=0> N_pool;
  int<lower=0> N_pred;
  real time_int[N_indiv + N_pool + N_pred];
  
  int<lower=0> k[N_indiv + N_pool];
  int<lower=0> y[N_indiv + N_pool];
  
  int<lower=0> num_group;
  int<lower=0> group[N_indiv];

  real<lower=0> ig_alpha;
  real<lower=0> ig_beta;
  int<lower=0> pool_size;
}

transformed data {
  real delta = 1e-9;
  int<lower=1> N = N_indiv + N_pool + N_pred;
  int<lower=0> N_pool_start = N_indiv + 1;
  int<lower=0> N_data = N_indiv + N_pool;
}

parameters {
  real<lower=0> elleq;
  real<lower=0> sigma;
  real<lower=0> sigma_mu;
  real mu;
  vector[N] eta;
  vector[num_group] mu_group;
}

transformed parameters {
  vector[N] z;
  vector[N] mu_vec;
  
  for (i in N_pool_start:N){
    mu_vec[i] = mu;
  }
  
  for (i in 1:N_indiv){
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
  y[1:N_indiv] ~ binomial(k[1:N_indiv], Phi(z[1:N_indiv])); //
   for (i in N_pool_start:N_data){
    y[i] ~ binomial(k[i], 1 - (1 - Phi(z[i]))^pool_size); 
  }
  mu_group ~ normal(mu, sigma_mu);
  sigma_mu ~ inv_gamma(1, 10);
}
