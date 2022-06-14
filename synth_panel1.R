library(tidyverse)
library(lubridate)
library(rstan)
library(mvtnorm)
library(shinystan)
# library(invgamma)

options(mc.cores = 4)

my_seed <- 370119 
set.seed(my_seed)

N <- 25 # num observed time points (weekly for 1 year)
t_range <- 1000 # num total days in date range
k = 15 # number of pools 
m = 5 # pool size

# True hyperparameters for exponentiated quadratic (squared exp) kernel
ell <- 100
sigma <- .5
mu01 <- 0.025

t <- seq(1, t_range - 7, by = 1)
dt <- plgp::distance(t)
eps <- sqrt(.Machine$double.eps) 

# Compute covariance
Sigma <- sigma^2 * exp(-dt / (2*ell^2)) + diag(eps, length(t)) 
# Simulate data
yprev <- rmvnorm(1, mean=rep(qnorm(mu01), length(t)), sigma = Sigma)

sim_df <- tibble(t = t, y = pnorm(drop(yprev)))

# Sample prevalence at N evenly spaced time points
observed <- seq(1, t_range - 7, by=(t_range %/% N))
obs_times <- t[c(observed)]
obs_prevs <- pnorm(yprev[1,c(observed)])

tnew <- seq(0, max(t), by=8)

# Plot prevalence curve and "observed" prevalence values
ggplot() + geom_line(aes(x=t, y=pnorm(yprev)), color='black') +
  geom_point(aes(x=obs_times, y=obs_prevs)) +
  theme_bw() + ylim(0,0.4)

pool_counts <- NULL
indiv_counts <- matrix(nrow = N, ncol = k * m)
indiv_obs_k <- NULL

for (i in 1:N){
  indiv_counts[i,] <- rbinom(1, n = m*k, prob = obs_prevs[i])
  pools <- matrix(indiv_counts[i,], ncol = m, byrow = T) %>% rowSums()
  indiv_obs_k[i] <- sum(indiv_counts[i,1:k])
  pool_counts[i] <- sum(pools > 0)
}

fit <- stan("stan_files/SIMpooled.stan",
            data = list(N1 = N,
                        N2 = length(tnew),
                        time_int = c(obs_times, tnew),
                        k = k,
                        m = m,
                        y = pool_counts,
                        ig_alpha = 7.3,
                        ig_beta = 675.2),
            control = list(adapt_delta = 0.9),
            iter = 3000, warmup = 1000,
            seed = my_seed)

# init = list(list(elleq=10), list(elleq=11), list(elleq=12), list(elleq=13)),

# fit2 <- stan("stan_files/SIMindiv.stan", 
#              data = list(N1 = N,
#                          N2 = length(tnew),
#                          time_int = c(obs_times, tnew),
#                          k = k,
#                          y = indiv_obs_k, 
#                          ig_alpha = 7.3,
#                          ig_beta = 675.2), 
#              control = list(adapt_delta = 0.9),
#              seed = my_seed)

# fit3 <- stan("stan_files/SIMindiv.stan", 
#              data = list(N1 = N,
#                          N2 = length(tnew),
#                          time_int = c(obs_times, tnew),
#                          k = k * m,
#                          y = rowSums(indiv_counts), 
#                          ig_alpha = 7.3,
#                          ig_beta = 675.2), 
#              control = list(adapt_delta = 0.9),
#              seed = my_seed)

# init = list(list(elleq=10), list(elleq=11), list(elleq=12), list(elleq=13)),

print(fit, pars=c("elleq", "sigma", "mu"), probs=c(0.025, 0.975))
print(fit2, pars=c("elleq", "sigma", "mu"), probs=c(0.025, 0.975))
print(fit3, pars=c("elleq", "sigma", "mu"), probs=c(0.025, 0.975))

zzz <- extract(fit)$z
zzz2 <- extract(fit2)$z
zzz3 <- extract(fit3)$z

out_pool <- tibble(t = c(obs_times, tnew),
                   mean = pnorm(apply(zzz[1:8000,], 2, mean)),
                   lower = pnorm(apply(zzz[1:8000,], 2, quantile, prob=0.025)),
                   uper = pnorm(apply(zzz[1:8000,], 2, quantile, prob=0.975)))

out_indiv <- tibble(t = c(obs_times, tnew),
                    mean = pnorm(apply(zzz2[1:4000,], 2, mean)),
                    lower = pnorm(apply(zzz2[1:4000,], 2, quantile, prob=0.025)),
                    uper = pnorm(apply(zzz2[1:4000,], 2, quantile, prob=0.975)))

out_indiv_km <- tibble(t = c(obs_times, tnew),
                    mean = pnorm(apply(zzz3[1:4000,], 2, mean)),
                    lower = pnorm(apply(zzz3[1:4000,], 2, quantile, prob=0.025)),
                    uper = pnorm(apply(zzz3[1:4000,], 2, quantile, prob=0.975)))

pool_obs <- tibble(t = obs_times,
                   freq = pool_counts / k)
indiv_obs <- tibble(t = obs_times,
                   freq = indiv_obs_k / k)
indiv_obs_km <- tibble(t = obs_times,
                    freq = rowSums(indiv_counts) / (k*m))


saveRDS(list(out_pool, out_indiv, pool_obs, indiv_obs, sim_df, out_indiv_km, indiv_obs_km), "SIMpanel1_370119")
my_list <- readRDS("SIMpanel1_370119")

ggplot() +
  geom_ribbon(data = my_list[[2]], aes(x=t, ymin=lower, ymax=uper), fill='blue', alpha=0.25) +
  geom_ribbon(data = my_list[[1]], aes(x=t, ymin=lower, ymax=uper), fill='red', alpha=0.5) +
  geom_line(aes(x=t, y=pnorm(yprev)), color='black') + theme_classic() +
  geom_line(data = my_list[[1]], aes(x=t, y=median), linetype=2) +
  geom_line(data = my_list[[2]], aes(x=t, y=median), linetype=3) +
  geom_point(data=my_list[[3]], aes(x=t, y=freq), shape=1) +
  geom_point(data=my_list[[4]], aes(x=t, y=freq), shape=4) +
  ylim(0,.5)
