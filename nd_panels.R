library(tidyverse)
library(lubridate)
library(rstan)

options(mc.cores = 4)

my_seed <- 5784
set.seed(my_seed)

dat <- read_csv("my_data/ND_data.csv") %>% filter(surveil_total > 0) 

t <- mdy(dat$date) %>% as.integer()
N <- length(dat$surveil_pos)

pool_sizes <- c(5, 10)

k1 <- dat$surveil_total %/% pool_sizes[1] 
m1 <- dat$surveil_total %% pool_sizes[1]

y1 <- matrix(nrow=N, ncol=2)
for (i in 1:N){
  final <- (k1[i] * pool_sizes[1])
  indiv_results <- c(rep(1, dat$surveil_pos[i]), 
                     rep(0, dat$surveil_total[i] -  dat$surveil_pos[i])) %>%
                   sample(., dat$surveil_total[i])
  if (k1[i] == 0){
    y1[i,1] <- 0
    y1[i,2] <- ifelse(sum(indiv_results) == 0, 0, 1)
  } else {
    pool_results <- matrix(indiv_results[1:final], nrow = k1[i], byrow = TRUE) %>% rowSums()
    y1[i,1] <- sum(ifelse(pool_results == 0, 0, 1))
    if (m1[i] == 0){
      y1[i,2] <- 0
    } else{
      y1[i,2] <- ifelse(sum(indiv_results[(final + 1):(final + m1[i])]) == 0, 0, 1)  
    }
  }
}

k2 <- dat$surveil_total %/% pool_sizes[2] 
m2 <- dat$surveil_total %% pool_sizes[2]

y2 <- matrix(nrow=N, ncol=2)
for (i in 1:N){
  final <- (k2[i] * pool_sizes[2])
  indiv_results <- c(rep(1, dat$surveil_pos[i]), 
                     rep(0, dat$surveil_total[i] -  dat$surveil_pos[i])) %>%
    sample(., dat$surveil_total[i])
  if (k2[i] == 0){
    y2[i,1] <- 0
    y2[i,2] <- ifelse(sum(indiv_results) == 0, 0, 1)
  } else {
    pool_results <- matrix(indiv_results[1:final], nrow = k2[i], byrow = TRUE) %>% rowSums()
    y2[i,1] <- sum(ifelse(pool_results == 0, 0, 1))
    if (m2[i] == 0){
      y2[i,2] <- 0
    } else{
      y2[i,2] <- ifelse(sum(indiv_results[(final + 1):(final + m2[i])]) == 0, 0, 1)  
    }
  }
}

pooled5 <- stan("stan_files/NDpooled.stan", 
                    data = list(N = N,
                                time_int = t,
                                k = k1,
                                m_all = pool_sizes[1],
                                m = m1,
                                y = y1,
                                ig_alpha = 11.8, 
                                ig_beta = 297.7),
                    iter = 1000, warmup = 500, chains=4, 
                    control = list(adapt_delta = 0.9, stepsize = 1, max_treedepth = 10),
                seed = my_seed) 

z_vals5 <- extract(pooled5)$z
out5 <- tibble(median = pnorm(apply(z_vals5[1:2000,], 2, median)),
               lower =  pnorm(apply(z_vals5[1:2000,], 2, quantile, prob = .025)),
               upper =  pnorm(apply(z_vals5[1:2000,], 2, quantile, prob = .975)), 
               x = t)

pooled10 <- stan("stan_files/NDpooled.stan", 
                    data = list(N = N,
                                time_int = t,
                                k = k2,
                                m_all = pool_sizes[2],
                                m = m2,
                                y = y2,
                                ig_alpha = 11.8, 
                                ig_beta = 297.7),
                    iter = 1000, warmup = 500, chains=4, 
                    init = list(list(elleq=20), list(elleq=21), list(elleq=22), list(elleq=23)),
                    control = list(adapt_delta = 0.9, max_treedepth = 10),
                 seed = my_seed) 

z_vals10 <- extract(pooled10)$z
out10 <- tibble(median = pnorm(apply(z_vals10[1:1500,], 2, median)),
               lower =  pnorm(apply(z_vals10[1:1500,], 2, quantile, prob = .025)),
               upper =  pnorm(apply(z_vals10[1:1500,], 2, quantile, prob = .975)), 
               x = t)

indivALL <- stan("stan_files/NDindiv.stan",
                   data=list(k = dat$surveil_total, 
                             N = length(t), 
                             y = dat$surveil_pos, 
                             time_int = t,
                             ig_alpha = 11.8, 
                             ig_beta = 297.7),
                   iter = 1000, warmup = 500, chains=4,
                   control = list(adapt_delta=0.8), 
                   seed = my_seed)

z_valsALL <- extract(indivALL)$z
outALL <- tibble(median = pnorm(apply(z_valsALL[1:2000,], 2, median)),
               lower =  pnorm(apply(z_valsALL[1:2000,], 2, quantile, prob = .025)),
               upper =  pnorm(apply(z_valsALL[1:2000,], 2, quantile, prob = .975)), 
               x = t)

y_indivLIMITED <- NULL
for (i in 1:N) {
  indiv_results <- c(rep(1, dat$surveil_pos[i]), rep(0, dat$surveil_total[i] - dat$surveil_pos[i]))
  y_indivLIMITED[i] <- sum(sample(indiv_results, k1[i] + ifelse(m1[i]==0, 0, 1), replace=FALSE))
}

indivLIMITED <- stan("stan_files/NDindiv.stan",
                 data=list(k = k1 + ifelse(m1 == 0, 0, 1), 
                           N = length(t), 
                           y = y_indivLIMITED, 
                           time_int = t,
                           ig_alpha = 11.8, 
                           ig_beta = 297.7),
                 iter = 1500, warmup = 750, chains=4,
                 control = list(adapt_delta=0.9), 
                 seed = my_seed)

print(pooled5, pars=c("elleq", "sigma", "mu"), probs=c(0.025, 0.975))
print(pooled10, pars=c("elleq", "sigma", "mu"), probs=c(0.025, 0.975))
print(indivALL, pars=c("elleq", "sigma", "mu"), probs=c(0.025, 0.975))
print(indivLIMITED, pars=c("elleq", "sigma", "mu"), probs=c(0.025, 0.975))


z_valsLIMITED <- extract(indivLIMITED)$z
outLIMITED <- tibble(median = pnorm(apply(z_valsLIMITED[1:3000,], 2, median)),
               lower =  pnorm(apply(z_valsLIMITED[1:3000,], 2, quantile, prob = .025)),
               upper =  pnorm(apply(z_valsLIMITED[1:3000,], 2, quantile, prob = .975)), 
               x = t)

# nd_out <- saveRDS(list(out5, out10, outALL, outLIMITED), file = "nd_out_5784")

ggplot() + 
  geom_point(aes(x=t, y=(dat$surveil_pos / dat$surveil_total))) + theme_bw() +
  geom_ribbon(data=out10, aes(x=x, ymin=lower, ymax=upper), linetype=1, 
              fill='blue', alpha=0.5) +
  geom_ribbon(data=out5, aes(x=x, ymin=lower, ymax=upper), linetype=1, 
              fill='red', alpha=.5) +
  geom_ribbon(data=outALL, aes(x=x, ymin=lower, ymax=upper), linetype=1, 
              color='black', alpha=0.) +
  geom_line(data=outALL, aes(x=x, y=median), linetype=1, color='black') +
  geom_line(data=out10, aes(x=x, y=median), linetype=1, color='red') + 
  geom_line(data=out5, aes(x=x, y=median), linetype=2, color='red') +
  ylim(0,0.06) + theme_bw()

