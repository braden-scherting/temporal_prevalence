library(tidyverse)
library(lubridate)
library(rstan)
library(mvtnorm)
library(shinystan)
library(invgamma)

options(mc.cores = 4)

my_seed <- 489713 
set.seed(my_seed)

N <- 150 # num observed time points (weekly for 1 year)
t_range <- 150 # num total days in date range
k = 50 # number of pools 
m = 10 # pool size

# True hyperparameters for exponentiated quadratic (squared exp) kernel
ell <- 20
sigma <- .25
mu01 <- 0.01

t <- 1:N
dt <- plgp::distance(t)
eps <- sqrt(.Machine$double.eps) 

# Compute covariance
Sigma <- sigma^2 * exp(-dt / (2*ell^2)) + diag(eps, length(t)) 
# Simulate data
yprev <- rmvnorm(1, mean=rep(qnorm(mu01), length(t)), sigma = Sigma)

sim_df <- tibble(t = t, y = pnorm(drop(yprev)))

# Sample prevalence at N evenly spaced time points
observed <- 1:N
obs_times <- t
obs_prevs <- pnorm(yprev)

# Plot prevalence curve and "observed" prevalence values
ggplot() + geom_line(aes(x=t, y=pnorm(yprev)), color='black') +
  geom_point(aes(x=obs_times, y=obs_prevs)) +
  theme_bw() + ylim(0,0.05)

pool_counts <- NULL
indiv_counts <- matrix(nrow = N, ncol = k * m)
indiv_obs_k <- NULL

for (i in 1:N){
  indiv_counts[i,] <- rbinom(1, n = m*k, prob = obs_prevs[i])
  pools <- matrix(indiv_counts[i,], ncol = m, byrow = T) %>% rowSums()
  indiv_obs_k[i] <- sum(sample(indiv_counts[i,], k))
  pool_counts[i] <- sum(pools > 0)
}


fit <- stan("stan_files/SIMpooled.stan", 
            data = list(N1 = N,
                        N2 = 0,
                        time_int = obs_times,
                        k = k,
                        m = m,
                        y = pool_counts, 
                        ig_alpha = 11.8, 
                        ig_beta = 297.7), 
            control = list(adapt_delta = 0.9),
            iter=2500, warmup = 500,
            seed = my_seed)

fit2 <- stan("stan_files/SIMindiv.stan", 
             data = list(N1 = N,
                         N2 = 0,
                         time_int = c(obs_times),
                         k = k,
                         y = indiv_obs_k, 
                         ig_alpha = 11.8, 
                         ig_beta = 297.7), 
             control = list(adapt_delta = 0.8),
             seed = my_seed)

fit3 <- stan("stan_files/SIMindiv.stan", 
             data = list(N1 = N,
                         N2 = 0,
                         time_int = c(obs_times),
                         k = k * m,
                         y = rowSums(indiv_counts), 
                         ig_alpha = 11.8, 
                         ig_beta = 297.7), 
             control = list(adapt_delta = 0.9),
             seed = my_seed)

print(fit, pars=c("elleq", "sigma", "mu"), probs=c(0.025, 0.975))
print(fit2, pars=c("elleq", "sigma", "mu"), probs=c(0.025, 0.975))
print(fit3, pars=c("elleq", "sigma", "mu"), probs=c(0.025, 0.975))

zzz <- extract(fit)$z
zzz2 <- extract(fit2)$z
zzz3 <- extract(fit3)$z

out_pool <- tibble(t = c(obs_times),
                   mean = pnorm(apply(zzz[1:8000,], 2, mean)),
                   lower = pnorm(apply(zzz[1:8000,], 2, quantile, prob=0.025)),
                   uper = pnorm(apply(zzz[1:8000,], 2, quantile, prob=0.975)))

out_indiv <- tibble(t = c(obs_times),
                    mean = pnorm(apply(zzz2[1:4000,], 2, mean)),
                    lower = pnorm(apply(zzz2[1:4000,], 2, quantile, prob=0.025)),
                    uper = pnorm(apply(zzz2[1:4000,], 2, quantile, prob=0.975)))

out_indiv_km <- tibble(t = c(obs_times),
                    mean = pnorm(apply(zzz3[1:4000,], 2, mean)),
                    lower = pnorm(apply(zzz3[1:4000,], 2, quantile, prob=0.025)),
                    uper = pnorm(apply(zzz3[1:4000,], 2, quantile, prob=0.975)))

pool_obs <- tibble(t = obs_times,
                   freq = pool_counts / k)
indiv_obs <- tibble(t = obs_times,
                    freq = indiv_obs_k / k)
indiv_obs_km <- tibble(t = obs_times,
                    freq = rowSums(indiv_counts) / (k*m))

saveRDS(list(out_pool, out_indiv, pool_obs, indiv_obs, sim_df, out_indiv_km,
             indiv_obs_km), "SIMpanel_489713")

my_list <- readRDS("SIMpanel_489713")

ggplot() +
  geom_ribbon(data = my_list[[2]], aes(x=t, ymin=lower, ymax=uper), fill='blue', alpha=0.25) +
  geom_ribbon(data = my_list[[1]], aes(x=t, ymin=lower, ymax=uper), fill='red', alpha=0.5) +
  geom_line(aes(x=t, y=pnorm(yprev)), color='black') +
  geom_line(data = my_list[[2]], aes(x=t, y=mean), color='black', linetype=2) +
  geom_line(data = my_list[[1]], aes(x=t, y=mean), color='black', linetype=3) +
  ylim(0,0.06) + theme_classic()
