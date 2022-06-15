library(tidyverse)
library(lubridate)
library(rstan)

options(mc.cores = 4)

my_seed <- 9052
set.seed(my_seed)

dat <- read_csv("my_data/batDataAggregated.csv")[23:42,]

t <- mdy(dat$date) %>% as.integer()
tnew <- seq(min(t-10), max(t + 10), by=4)
n_pred <- length(tnew)
N <- length(dat$n_pos)

pool_sizes <- c(3, 5)
k1 <- dat$n_tests %/% pool_sizes[1] + ifelse(dat$n_tests %% pool_sizes[1] == 0, 0, 1)
m1 <- matrix(rep(0, N * max(k1)), nrow=N)
for (i in 1:N){
  m1[i,1:(k1[i] - 1)] <- pool_sizes[1]
  if (dat$n_tests[i] %% pool_sizes[1] == 0){
    m1[i,k1[i]] <- pool_sizes[1]
  } else {
    m1[i, k1[i]] <- dat$n_tests[i] %% pool_sizes[1]  
  }
}

y1 <- matrix(rep(0, N * max(k1)), nrow=N)
for (i in 1:N){
  indiv_results <- c(rep(1, dat$n_pos[i]), rep(0, dat$n_tests[i] - dat$n_pos[i]))
  pool_results <- c(sample(indiv_results, dat$n_tests[i]), 
                    rep(0, k1[i]*pool_sizes[1] - length(indiv_results))) %>%
    matrix(., nrow=k1[i], byrow=TRUE) %>% rowSums()
  y1[i,1:k1[i]] <- ifelse(pool_results == 0, 0, 1)
}

k2 <- dat$n_tests %/% pool_sizes[2] + ifelse(dat$n_tests %% pool_sizes[2] == 0, 0, 1)
m2 <- matrix(rep(0, N * max(k2)), nrow=N)
for (i in 1:N){
  m2[i,1:(k2[i] - 1)] <- pool_sizes[2]
  if (dat$n_tests[i] %% pool_sizes[2] == 0){
    m2[i,k2[i]] <- pool_sizes[2]
  } else {
    m2[i, k2[i]] <- dat$n_tests[i] %% pool_sizes[2]  
  }
}

y2 <- matrix(rep(0, N * max(k2)), nrow=N)
for (i in 1:N){
  indiv_results <- c(rep(1, dat$n_pos[i]), rep(0, dat$n_tests[i] - dat$n_pos[i]))
  pool_results <- c(sample(indiv_results, dat$n_tests[i]), 
                    rep(0, k2[i]*pool_sizes[2] - length(indiv_results))) %>%
    matrix(., nrow=k2[i], byrow=TRUE) %>% rowSums()
  y2[i,1:k2[i]] <- ifelse(pool_results == 0, 0, 1)
}

pooled3 <- stan("stan_files/BATpooled.stan",
            data = list(N1 = N,
                        N2 = length(tnew),
                        time_int = c(t, tnew),
                        max_k = max(k1),
                        k = k1,
                        m = m1,
                        y = y1,
                        ig_alpha = 7.3,
                        ig_beta = 675.2),
            iter = 1200, chains=4,
            seed = my_seed * pool_sizes[1])
z_vals3 <- extract(pooled3)$z
out3 <- tibble(mean = pnorm(apply(z_vals3[1:2400,], 2, mean)),
               lower =  pnorm(apply(z_vals3[1:2400,], 2, quantile, prob = .025)),
               upper =  pnorm(apply(z_vals3[1:2400,], 2, quantile, prob = .975)), 
               x = c(t,tnew))

pooled5 <- stan("stan_files/BATpooled.stan",
                data = list(N1 = N,
                            N2 = length(tnew),
                            time_int = c(t, tnew),
                            max_k = max(k2),
                            k = k2,
                            m = m2,
                            y = y2,
                            ig_alpha = 7.3,
                            ig_beta = 675.2),
                iter = 1200, chains=4,
                seed = my_seed * pool_sizes[2])
z_vals5 <- extract(pooled5)$z
out5 <- tibble(mean = pnorm(apply(z_vals5[1:2400,], 2, mean)),
               lower =  pnorm(apply(z_vals5[1:2400,], 2, quantile, prob = .025)),
               upper =  pnorm(apply(z_vals5[1:2400,], 2, quantile, prob = .975)), 
               x = c(t,tnew))


indivALL <- stan("stan_files/BATindiv.stan",
                     data = list(N1 = N,
                                 N2 = length(tnew),
                                 time_int = c(t, tnew),
                                 k = dat$n_tests,
                                 y = dat$n_pos,
                                 ig_alpha = 7.3,
                                 ig_beta = 675.2),
                     iter = 1200, chains=4,
                     seed = my_seed)
z_valsALL <- extract(indivALL)$z
outALL <- tibble(mean = pnorm(apply(z_valsALL[1:2400,], 2, mean)),
                     lower =  pnorm(apply(z_valsALL[1:2400,], 2, quantile, prob = .025)),
                     upper =  pnorm(apply(z_valsALL[1:2400,], 2, quantile, prob = .975)), 
                     x = c(t, tnew))


y_indivLIMITED <- NULL
for (i in 1:N) {
  indiv_results <- c(rep(1, dat$n_pos[i]), rep(0, dat$n_tests[i] - dat$n_pos[i]))
  y_indivLIMITED[i] <- sum(sample(indiv_results, k1[i], replace=FALSE))
}

indivLIMITED <- stan("stan_files/BATindiv.stan",
                     data = list(N1 = N,
                                 N2 = length(tnew),
                                 time_int = c(t, tnew),
                                 k = k1,
                                 y = y_indivLIMITED,
                                 ig_alpha = 7.3,
                                 ig_beta = 675.2),
                     iter = 1500, chains=4,
                     seed = my_seed)
z_valsLIMITED <- extract(indivLIMITED)$z
outLIMITED <- tibble(mean = pnorm(apply(z_valsLIMITED[1:2400,], 2, mean)),
               lower =  pnorm(apply(z_valsLIMITED[1:2400,], 2, quantile, prob = .025)),
               upper =  pnorm(apply(z_valsLIMITED[1:2400,], 2, quantile, prob = .975)), 
               x = c(t,tnew))

print(pooled3, pars=c("elleq", "sigma", "mu"), probs=c(0.025, 0.975))
print(pooled5, pars=c("elleq", "sigma", "mu"), probs=c(0.025, 0.975))
print(indivALL, pars=c("elleq", "sigma", "mu"), probs=c(0.025, 0.975))
print(indivLIMITED, pars=c("elleq", "sigma", "mu"), probs=c(0.025, 0.975))

bat_out <- saveRDS(list(out3, out5, outALL, outLIMITED), file = "bat_out_9052")
bat_out <- readRDS("bat_out_9052")
ggplot() + 
  geom_point(aes(x=as.Date(t,origin), y=(dat$n_pos / dat$n_tests))) + theme_bw() +
  geom_ribbon(data=outLIMITED, aes(x=as.Date(x,origin), ymin=lower, ymax=upper), linetype=1, 
              fill='yellow', alpha=0.5, color='black') +
  geom_ribbon(data=out5, aes(x=as.Date(x,origin), ymin=lower, ymax=upper), linetype=1, 
              fill='blue', alpha=0.5) +
  geom_ribbon(data=out3, aes(x=as.Date(x,origin), ymin=lower, ymax=upper), linetype=1, 
              fill='red', alpha=.5) +
  geom_ribbon(data=outALL, aes(x=as.Date(x,origin), ymin=lower, ymax=upper), linetype=1, 
              color='black', alpha=0.) +
  geom_line(data=outALL, aes(x=as.Date(x,origin), y=mean), linetype=1, color='black') +
  geom_line(data=out3, aes(x=as.Date(x,origin), y=mean), linetype=1, color='red') + 
  geom_line(data=out5, aes(x=as.Date(x,origin), y=mean), linetype=2, color='red') +
  geom_line(data=outLIMITED, aes(x=as.Date(x,origin), y=mean), linetype=2, color='yellow') +
  ylim(0,0.7) + theme_bw()
