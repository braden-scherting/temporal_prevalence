library(tidyverse)
library(lubridate)
library(rstan)
library(mvtnorm)
library(shinystan)
library(tictoc)

options(mc.cores = 4)

my_seed <- 03262022 
set.seed(my_seed)

N <- 25 # num observed time points (weekly for 1 year)
t_range <- 1000 # num total days in date range
k = 15 # number of pools 
m = 5 # pool size

# True hyperparameters for exponentiated quadratic (squared exp) kernel
ell <- 100
sigma <- .5
mu01 <- qnorm(0.025)

num_groups <- 5
sd_mu <- .2
mu <- rnorm(num_groups, mu01, sd_mu)

t <- seq(1, t_range - 7, by = 1)
dt <- plgp::distance(t)
eps <- sqrt(.Machine$double.eps) 

# Compute covariance
Sigma <- sigma^2 * exp(-dt / (2*ell^2)) + diag(eps, length(t)) 
# Simulate data
eta <- rmvnorm(1, mean=rep(0, length(t)), sigma = Sigma)
yprev <- eta + mu01
# yprev <- rmvnorm(1, mean=rep(qnorm(mu01), length(t)), sigma = Sigma)

sim_df <- tibble(time = rep(t, num_groups + 1), 
                 y = c(pnorm(drop(yprev)), pnorm(drop(eta) + mu[1]) ,
                       pnorm(drop(eta) + mu[2]), pnorm(drop(eta) + mu[3]),
                       pnorm(drop(eta) + mu[4]),
                       pnorm(drop(eta) + mu[5])
                       #pnorm(drop(yprev) + mu[6])
                       #, pnorm(drop(yprev) + mu[7]),
                       # pnorm(drop(yprev) + mu[8]), pnorm(drop(yprev) + mu[9]),
                       #  pnorm(drop(yprev) + mu[10])
                 ),
                 group = rep(c('Overall', 1:num_groups), each = length(t)))



# Sample prevalence at N evenly spaced time points
observed <- seq(1, t_range - 7, by=(t_range %/% N))
obs_times <- t[c(observed)]
obs_prevs <- pnorm(yprev[1,c(observed)])

group_prevs <- cbind(pnorm((eta + mu[1])[1,c(observed)]),
                     pnorm((eta + mu[2])[1,c(observed)]),
                     pnorm((eta + mu[3])[1,c(observed)]),
                     pnorm((eta + mu[4])[1,c(observed)]),
                     pnorm((eta + mu[5])[1,c(observed)])
                     #pnorm((yprev + mu[6])[1,c(observed)])#,
                     #  pnorm((yprev + mu[7])[1,c(observed)]),
                     # pnorm((yprev + mu[8])[1,c(observed)]),
                     #pnorm((yprev + mu[9])[1,c(observed)]),
                     #pnorm((yprev + mu[10])[1,c(observed)])
)

tnew <- seq(0, max(t), by=8)

# Plot prevalence curve and "observed" prevalence values
ggplot() + geom_line(aes(x=t, y=pnorm(yprev)), color='black') +
  geom_point(aes(x=obs_times, y=obs_prevs)) +
  theme_bw() + ylim(0,0.5)

# Plot prevalence curve and "observed" prevalence values
sim_df %>% ggplot(aes(x=time, y= y, color = group)) + 
  geom_line() +
  geom_line(aes(x=time, y=y), color='black', inherit.aes = F, data = sim_df %>% filter(group == "Overall") )+
  # geom_point(aes(x=obs_times, y=obs_prevs), inherit.aes = F, data = sim_df %>% filter(group == "Overall")) +
  theme_bw() +  theme(legend.position="none") +
  ylab('Prevalence') + ggtitle("Latent Prevalence") + 
  ylim(0,0.5)

pool_counts <- NULL
indiv_counts <- matrix(nrow = N, ncol = k * m)
indiv_obs_k <- NULL

for (i in 1:N){
  group_ids <- sample(num_groups, m*k, replace = T)
  indiv_counts[i,] <- rbinom(1, n = m*k, prob = group_prevs[i, group_ids])
  pools <- matrix(indiv_counts[i,], ncol = m, byrow = T) %>% rowSums()
  pool_counts[i] <- sum(pools > 0)
}

num_reps <- 100
indiv_counts <- matrix(nrow = N, ncol = num_groups * num_reps)
indiv_obs_k <- NULL

for (i in 1:N){
  indiv_counts[i,] <- rbinom(1, n = num_groups * num_reps, prob = group_prevs[i, rep(1:num_groups, each = num_reps)])
  indiv_obs_k[i] <- sum(indiv_counts[i,])
}

counts <- tibble(counts = as.numeric(t(indiv_counts)), 
                 time = rep(1:N, each = num_groups * num_reps), 
                 group = rep(rep(1:num_groups, each = num_reps), N)) %>%
  group_by(time,group) %>%
  summarize(counts = sum(counts), .groups = 'drop') %>%
  select(counts) %>% pull()
# 
# ################################################
# #  POOLED DATA
# ###### random sample from each species
# ################################################
# 
# tic()
# fit <- stan("stan_files/SIMpooled.stan",
#             data = list(N1 = N,
#                         N2 = length(tnew),
#                         time_int = c(obs_times, tnew),
#                         k = k,
#                         m = m,
#                         y = pool_counts,
#                         ig_alpha = 7.3,
#                         ig_beta = 675.2),
#             control = list(adapt_delta = 0.9),
#             iter = 3000, warmup = 1000,
#             seed = my_seed)
# toc()
# print(fit, pars=c("elleq", "sigma", "mu"), probs=c(0.025, 0.975))
# zzz <- extract(fit)$z
# 
# out_pool <- tibble(t = c(obs_times, tnew),
#                    median = pnorm(apply(zzz[1:8000,], 2, median)),
#                    lower = pnorm(apply(zzz[1:8000,], 2, quantile, prob=0.025)),
#                    uper = pnorm(apply(zzz[1:8000,], 2, quantile, prob=0.975)))
# pool_obs <- tibble(t = obs_times,
#                    freq = pool_counts / k)
# 
# ggplot() +
#   geom_ribbon(data = out_pool, aes(x=t, ymin=lower, ymax=uper), fill='red', alpha=0.5) +
#   geom_line(aes(x=t, y=pnorm(yprev)), color='black') + theme_bw() +
#   geom_line(data = out_pool, aes(x=t, y=median), linetype=2) +
#   ylim(0,.5) + ylab("Prevalence") +
#   ggtitle("Pooled Prevalence Estimates")
# 
# 
# ################################################
# #  INDIVIDUAL DATA
# ###### random sample from each species
# ###### analyze with binomial and no group label
# ################################################
# 
# tic()
# fit2 <- stan("stan_files/SIMindiv.stan",
#              data = list(N1 = N,
#                          N2 = length(tnew),
#                          time_int = c(obs_times, tnew),
#                          k = num_groups * num_reps,
#                          y = indiv_obs_k,
#                          ig_alpha = 7.3,
#                          ig_beta = 675.2),
#              control = list(adapt_delta = 0.9),
#              seed = my_seed)
# toc()
# 
# print(fit2, pars=c("elleq", "sigma", "mu"), probs=c(0.025, 0.975))
# zzz2 <- extract(fit2)$z
# out_indiv <- tibble(t = c(obs_times, tnew),
#                        median = pnorm(apply(zzz2[1:4000,], 2, median)),
#                        lower = pnorm(apply(zzz2[1:4000,], 2, quantile, prob=0.025)),
#                        uper = pnorm(apply(zzz2[1:4000,], 2, quantile, prob=0.975)))
# 
# indiv_obs <- tibble(t = obs_times, freq = indiv_obs_k / (num_groups * num_reps))
# ggplot() +
#   geom_ribbon(data = out_indiv, aes(x=t, ymin=lower, ymax=uper), fill='blue', alpha=0.25) +
#   geom_line(aes(x=t, y=pnorm(yprev)), color='black') + theme_classic() +
#   geom_line(data = out_indiv, aes(x=t, y=median), linetype=3) +
#   geom_point(data=indiv_obs, aes(x=t, y=freq), shape=4) +
#   ylim(0,.5)
# 
# # ################################################
# # #  INDIVIDUAL DATA
# # ###### random sample from each species
# # ###### analyze with Bernoulli and no group label
# # ################################################
# # 
# # tic()
# # fit3 <- stan("stan_files/SIMindiv.stan",
# #              data = list(N1 = N * num_groups * num_reps,
# #                          N2 = length(tnew),
# #                          time_int = c(rep(obs_times, each = num_groups * num_reps), tnew),
# #                          k = 1,
# #                          y = as.numeric(t(indiv_counts)),
# #                          ig_alpha = 7.3,
# #                          ig_beta = 675.2),
# #              control = list(adapt_delta = 0.9),
# #              seed = my_seed)
# # toc()
# # 
# # print(fit3, pars=c("elleq", "sigma", "mu"), probs=c(0.025, 0.975))
# # zzz3 <- extract(fit3)$z
# # out_indiv <- tibble(t = c(rep(obs_times, each = num_groups * num_reps), tnew),
# #                     median = pnorm(apply(zzz3[1:4000,], 2, median)),
# #                     lower = pnorm(apply(zzz3[1:4000,], 2, quantile, prob=0.025)),
# #                     uper = pnorm(apply(zzz3[1:4000,], 2, quantile, prob=0.975)))
# # 
# # indiv_obs <- tibble(t = obs_times,
# #                     freq = indiv_obs_k / (num_groups * num_reps))
# # ggplot() +
# #   geom_ribbon(data = out_indiv, aes(x=t, ymin=lower, ymax=uper), fill='blue', alpha=0.25) +
# #   geom_line(aes(x=t, y=pnorm(yprev)), color='black') + theme_classic() +
# #   geom_line(data = out_indiv, aes(x=t, y=median), linetype=3) +
# #   geom_point(data=indiv_obs, aes(x=t, y=freq), shape=4) +
# #   ylim(0,.5)
# 
# ################################################
# #  INDIVIDUAL DATA
# ###### random sample from each species
# ###### analyze with grouped Bernoulli but no group label
# ################################################
# 
# tic()
# fit3 <- stan("stan_files/SIMindiv.stan",
#              data = list(N1 = N * num_groups,
#                          N2 = length(tnew),
#                          time_int = c(rep(obs_times, each = num_groups), tnew),
#                          k = num_reps,
#                          y = counts,
#                          ig_alpha = 7.3,
#                          ig_beta = 675.2),
#              control = list(adapt_delta = 0.9),
#              seed = my_seed)
# toc()
# 
# print(fit3, pars=c("elleq", "sigma", "mu"), probs=c(0.025, 0.975))
# zzz3 <- extract(fit3)$z
# out_indiv <- tibble(t = c(rep(obs_times, each = num_groups ), tnew),
#                     median = pnorm(apply(zzz3[1:4000,], 2, median)),
#                     lower = pnorm(apply(zzz3[1:4000,], 2, quantile, prob=0.025)),
#                     uper = pnorm(apply(zzz3[1:4000,], 2, quantile, prob=0.975)))
# 
# indiv_obs <- tibble(t = obs_times,
#                     freq = indiv_obs_k / (num_groups * num_reps))
# ggplot() +
#   geom_ribbon(data = out_indiv, aes(x=t, ymin=lower, ymax=uper), fill='blue', alpha=0.25) +
#   geom_line(aes(x=t, y=pnorm(yprev)), color='black') + theme_classic() +
#   geom_line(data = out_indiv, aes(x=t, y=median), linetype=3) +
#   geom_point(data=indiv_obs, aes(x=t, y=freq), shape=4) +
#   ylim(0,.5)
# 
# ################################################
# #  INDIVIDUAL DATA
# ###### random sample from each species
# ###### analyze hierarchical
# ################################################
# 
# tic()
# fit4 <- stan("stan_files/SIMhier.stan",
#              data = list(N1 = N * num_groups,
#                          N2 = length(tnew),
#                          time_int = c(rep(obs_times, each = num_groups), tnew),
#                          k = num_reps,
#                          y = counts,
#                          ig_alpha = 7.3,
#                          ig_beta = 675.2,
#                          num_group = num_groups,
#                          group = rep(1:num_groups, N)),
#              control = list(adapt_delta = 0.9),
#              seed = my_seed)
# toc()
# 
# print(fit4, pars=c("elleq", "sigma", "mu", "mu_group"), probs=c(0.025, 0.975))
# zzz4 <- extract(fit4)$z
# out_indiv <- tibble(t = c(rep(obs_times, each = num_groups), tnew),
#                     group_ids = c(rep(1:4, length(obs_times)), rep('mu', length(tnew))),
#                     median = pnorm(apply(zzz4[1:4000,], 2, median)),
#                     lower = pnorm(apply(zzz4[1:4000,], 2, quantile, prob=0.025)),
#                     uper = pnorm(apply(zzz4[1:4000,], 2, quantile, prob=0.975)))
# 
# indiv_obs <- tibble(counts = as.numeric(t(indiv_counts)) / num_reps, 
#                     time = rep(obs_times, each = num_groups * num_reps), 
#                     group_ids = rep(rep(1:num_groups, each = num_reps), N)) %>%
#   group_by(time,group_ids) %>%
#   summarize(counts = sum(counts), .groups = 'drop')
# 
# mu_data <- tibble(x=rep(t,5), 
#                   y=c(pnorm(drop(yprev)), pnorm(drop(yprev) + mu[1]), pnorm(drop(yprev) + mu[2]),
#                       pnorm(drop(yprev) + mu[3]), pnorm(drop(yprev) + mu[4])), 
#                   group_ids = rep(c('mu', 1:4), each = length(t)))
# 
# mu_data <- tibble(x=t,y=c(pnorm(drop(yprev))),group_ids = rep(c('mu'), each = length(t)))
# mu_data1 <- tibble(x=t,y=c(pnorm(drop(yprev)+ mu[1])),group_ids = rep(c('1'), each = length(t)))
# mu_data2 <- tibble(x=t,y=c(pnorm(drop(yprev)+ mu[2])),group_ids = rep(c('2'), each = length(t)))
# mu_data3 <- tibble(x=t,y=c(pnorm(drop(yprev)+ mu[3])),group_ids = rep(c('3'), each = length(t)))
# mu_data4 <- tibble(x=t,y=c(pnorm(drop(yprev)+ mu[4])),group_ids = rep(c('4'), each = length(t)))
# 
# 
# out_indiv %>% ggplot() +
#   geom_ribbon(aes(x=t, ymin=lower, ymax=uper), fill='blue', alpha=0.25) +
#   geom_line(aes(x=t, y=y), color='black', data = mu_data) + theme_bw() +
#   geom_line(aes(x=t, y=y), color='black', data = mu_data1) +
#   geom_line(aes(x=t, y=y), color='black', data = mu_data2) +
#   geom_line(aes(x=t, y=y), color='black', data = mu_data3) +
#   geom_line(aes(x=t, y=y), color='black', data = mu_data4) +
#   facet_wrap(~group_ids) + geom_line(data = out_indiv, aes(x=t, y=median), linetype=3) +
#   geom_point(data=indiv_obs, aes(x=time, y=counts), shape=4) +
#   ylim(0,.5) 


################################################
#  INDIVIDUAL + POOLed DATA
###### random sample from each species
###### analyze hierarchical
################################################

# tic()
# fit5 <- stan("stan_files/SIM_dataintegrate.stan",
#              data = list(N_indiv = N * num_groups,
#                          N_pool = N,
#                          N_pred = length(tnew),
#                          time_int = c(rep(obs_times, each = num_groups),obs_times, tnew),
#                          k = c(rep(num_reps, N * num_groups), rep(k, N)),
#                          y = c(counts, pool_counts),
#                          ig_alpha = 7.3,
#                          ig_beta = 675.2,
#                          num_group = num_groups,
#                          group = rep(1:num_groups, N),
#                          pool_size = m),
#              control = list(adapt_delta = 0.9),
#              seed = my_seed)
# toc()
# 
# save(fit5, file = "integrate.RData")
load("integrate.RData")
# launch_shinystan(fit5)
print(fit5, pars=c("elleq", "sigma", "mu", "mu_group"), probs=c(0.025, 0.975))
zzz5 <- extract(fit5)$z
out_indiv <- tibble(t = c(rep(obs_times, each = num_groups), obs_times, tnew),
                    group_ids = c(rep(paste('group', 1:5), length(obs_times)), 
                                  rep('overall', length(obs_times)), 
                                  rep('overall', length(tnew))),
                    mean = pnorm(apply(zzz5[1:4000,], 2, mean)),
                    lower = pnorm(apply(zzz5[1:4000,], 2, quantile, prob=0.025)),
                    uper = pnorm(apply(zzz5[1:4000,], 2, quantile, prob=0.975)))

indiv_obs <- tibble(counts = as.numeric(t(indiv_counts)) / num_reps, 
                    time = rep(obs_times, each = num_groups * num_reps), 
                    group_ids = rep(rep(paste('group', 1:5), each = num_reps), N)) %>%
  group_by(time,group_ids) %>%
  summarize(counts = sum(counts), .groups = 'drop')

mu_data <- tibble(x=rep(t,6), 
                  y=c(pnorm(drop(yprev)), pnorm(drop(yprev) + mu[1]), pnorm(drop(yprev) + mu[2]),
                      pnorm(drop(yprev) + mu[3]), pnorm(drop(yprev) + mu[4]), pnorm(drop(yprev) + mu[5])), 
                  group_ids = rep(c('overall', paste('group', 1:5)), each = length(t)))

mu_data <- tibble(x=t,y=c(pnorm(drop(eta) + mu01)),group_ids = rep(c('overall'), each = length(t)))
mu_data1 <- tibble(x=t,y=c(pnorm(drop(eta) + mu[1])),group_ids = rep(c('group 1'), each = length(t)))
mu_data2 <- tibble(x=t,y=c(pnorm(drop(eta) + mu[2])),group_ids = rep(c('group 2'), each = length(t)))
mu_data3 <- tibble(x=t,y=c(pnorm(drop(eta) + mu[3])),group_ids = rep(c('group 3'), each = length(t)))
mu_data4 <- tibble(x=t,y=c(pnorm(drop(eta) + mu[4])),group_ids = rep(c('group 4'), each = length(t)))
mu_data5 <- tibble(x=t,y=c(pnorm(drop(eta) + mu[5])),group_ids = rep(c('group 5'), each = length(t)))


rug_data <- tibble(t = rep(obs_times, 6), group_ids = rep(c(paste('group', 1:5), 'overall'), each = length(obs_times))) 

pdf(file = "JSSAM/Figure5.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 6)
out_indiv %>% ggplot() +
  geom_ribbon(aes(x=t, ymin=lower, ymax=uper), fill='blue', alpha=0.25) +
  geom_line(aes(x=t, y=y), color='black', data = mu_data, linetype=3) + theme_bw() +
  geom_line(aes(x=t, y=y), color='black', data = mu_data1, linetype=3) +
  geom_line(aes(x=t, y=y), color='black', data = mu_data2, linetype=3) +
  geom_line(aes(x=t, y=y), color='black', data = mu_data3, linetype=3) +
  geom_line(aes(x=t, y=y), color='black', data = mu_data4, linetype=3) +
  geom_line(aes(x=t, y=y), color='black', data = mu_data5, linetype=3) +
  facet_wrap(~group_ids) + geom_line(data = out_indiv, aes(x=t, y=mean)) +
  geom_point(data=indiv_obs, aes(x=time, y=counts), shape=4) +
  ylim(0,.3) + geom_rug(aes(x = t), data = rug_data) + ylab("Prevalence") +
  xlab('Days') + ggtitle('Simulation study 3')

dev.off()
# plot(fit5, pars = c('mu', 'mu_group')) + 
#   geom_point(y =1 ,x =mu[5] , size = 3, pch = 23, fill = 'blue') +
#   geom_point(y =2 ,x =mu[4] , size = 3, pch = 23, fill = 'blue') +
#   geom_point(y =3 ,x =mu[3] , size = 3, pch = 23, fill = 'blue') +
#   geom_point(y =4 ,x =mu[2] , size = 3, pch = 23, fill = 'blue') +
#   geom_point(y =5 ,x =mu[1] , size = 3, pch = 23, fill = 'blue') +
#   geom_point(y =6 ,x =mu01 , size = 3, pch = 23, fill = 'blue') 
