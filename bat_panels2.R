library(tidyverse)
library(lubridate)
library(rstan)
library(tictoc)

options(mc.cores = 4)

my_seed <- 03302022
set.seed(my_seed)

bat_dat <- read.csv("my_data/batCoV.csv") %>% 
  mutate(Sampling_Date = mdy(Sampling_Date), positive = PCR_Result == 'Positive') %>%
  filter(Taxa %in% c('Bats','bats'), 
         year(Sampling_Date) < 2015,
         year(Sampling_Date) > 2010)


bat_dat <- bat_dat %>% group_by(Animal_ID) %>% sample_n(1) %>% ungroup()

bat_dat %>% group_by(Sampling_Date) %>% summarize(prev = mean(positive), n = n()) %>% 
  ggplot(aes(y = prev, x = Sampling_Date)) + geom_point() + geom_line()



species_interest <- bat_dat %>% group_by(Species) %>%
  summarize(num_dates = length(unique(Sampling_Date))) %>%
  arrange(desc(num_dates)) %>% slice(1:5) %>% select(Species) %>% pull()

bat_dat <- bat_dat %>% mutate(pool = !(Species %in% species_interest),
                              Species = case_when(pool ~ 'pool', TRUE ~ Species )) %>%
  select(Species, Sampling_Date, positive, pool)


# aggregated_data <- bat_dat %>% group_by(Sampling_Date, Species) %>%
#   summarise(pos = sum(positive), n = n(), .groups = 'drop') %>%
#   mutate(group = case_when(
#     Species == species_interest[1] ~ 1,
#     Species == species_interest[2] ~ 2,
#     Species == species_interest[3] ~ 3,
#     Species == species_interest[4] ~ 4,
#     Species == species_interest[5] ~ 5,
#     TRUE ~ 6),
#     julian_date = as.numeric(Sampling_Date))

indiv <- bat_dat %>% filter(!pool) %>% group_by(Sampling_Date, Species) %>%
  summarise(pos = sum(positive), n = n(), .groups = 'drop') %>%
  mutate(group = case_when(
    Species == species_interest[1] ~ 1,
    Species == species_interest[2] ~ 2,
    Species == species_interest[3] ~ 3,
    Species == species_interest[4] ~ 4,
    Species == species_interest[5] ~ 5,
    TRUE ~ 6),
    julian_date = as.numeric(Sampling_Date))

pool <- bat_dat %>% filter(pool) %>% group_by(Sampling_Date) %>%
  mutate(pool_numb = sample(ceiling(n()/3), n(), replace = T)) %>% 
  group_by(Sampling_Date, Species, pool_numb) %>% summarise(pos = sum(positive), n = n(), .groups = 'drop')  %>%
  mutate(group = case_when(
    Species == species_interest[1] ~ 1,
    Species == species_interest[2] ~ 2,
    Species == species_interest[3] ~ 3,
    Species == species_interest[4] ~ 4,
    Species == species_interest[5] ~ 5,
    TRUE ~ 6),
    julian_date = as.numeric(Sampling_Date))


t <- c(pool$julian_date, indiv$julian_date)
tnew <- seq(min(t-10), max(t + 10), by=30)
N_pred <- length(tnew)

N_indiv <- indiv %>% nrow()
N_pool <- pool %>% nrow()

time_int <- c(indiv$julian_date, pool$julian_date, tnew)
k <- c(indiv$n, rep(1,N_pool))
y <- c(indiv$pos, as.numeric(pool$pos > 0))
num_group <- length(species_interest)
group <- indiv$group
pool_size <- c(rep(1,N_indiv), pool$n)

tic()
bat_fit <- stan("stan_files/BAT_dataintegrate.stan",
             data = list(N_indiv = N_indiv,
                         N_pool = N_pool,
                         N_pred = N_pred,
                         time_int = time_int,
                         k = k,
                         y = y,
                         ig_alpha = 7.3,
                         ig_beta = 675.2,
                         num_group = num_group,
                         group = group,
                         pool_size = pool_size),
             control = list(adapt_delta = 0.9),
             seed = my_seed)
toc()


save(bat_fit, file = "bat_integrate.RData")

#launch_shinystan(bat_fit)
print(bat_fit, pars=c("elleq", "sigma", "mu", "mu_group"), probs=c(0.025, 0.975))
zzz <- extract(bat_fit)$z
# out_indiv <- tibble(t = as_date(time_int),
#                     group_ids = c(indiv$Species, 
#                                   rep('mu', N_pool), 
#                                   rep('mu', N_pred)),
#                     median = pnorm(apply(zzz, 2, median)),
#                     lower = pnorm(apply(zzz, 2, quantile, prob=0.025)),
#                     uper = pnorm(apply(zzz, 2, quantile, prob=0.975)))

out_indiv <- tibble(t = as_date(time_int)[(N_indiv+1):(N_indiv+N_pool +N_pred)],
                    group_ids = c(rep('overall', N_pool), 
                                  rep('overall', N_pred)),
                    median = pnorm(apply(zzz[,(N_indiv+1):(N_indiv+N_pool +N_pred)], 2, median)),
                    lower = pnorm(apply(zzz[,(N_indiv+1):(N_indiv+N_pool +N_pred)], 2, quantile, prob=0.025)),
                    uper = pnorm(apply(zzz[,(N_indiv+1):(N_indiv+N_pool +N_pred)], 2, quantile, prob=0.975)))


eta <- zzz[,(N_indiv + N_pool + 1):(N_indiv + N_pool + N_pred)] - matrix(extract(bat_fit)$mu, nrow = 4000, ncol = N_pred)

zzz1 <- eta + matrix(extract(bat_fit)$mu_group[,1], nrow = 4000, ncol = N_pred)
zzz2 <- eta + matrix(extract(bat_fit)$mu_group[,2], nrow = 4000, ncol = N_pred)
zzz3 <- eta + matrix(extract(bat_fit)$mu_group[,3], nrow = 4000, ncol = N_pred)
zzz4 <- eta + matrix(extract(bat_fit)$mu_group[,4], nrow = 4000, ncol = N_pred)
zzz5 <- eta + matrix(extract(bat_fit)$mu_group[,5], nrow = 4000, ncol = N_pred)

preds_grp1 <- tibble(t = rep(as_date(tnew), 5),
                     median = c(pnorm(apply(zzz1, 2, median)),
                                pnorm(apply(zzz2, 2, median)),
                                pnorm(apply(zzz3, 2, median)),
                                pnorm(apply(zzz4, 2, median)),
                                pnorm(apply(zzz5, 2, median))),
                     lower = c(pnorm(apply(zzz1, 2, quantile, prob=0.025)),
                               pnorm(apply(zzz2, 2, quantile, prob=0.025)),
                               pnorm(apply(zzz3, 2, quantile, prob=0.025)),
                               pnorm(apply(zzz4, 2, quantile, prob=0.025)),
                               pnorm(apply(zzz5, 2, quantile, prob=0.025))),
                     uper = c(pnorm(apply(zzz1, 2, quantile, prob=0.975)),
                              pnorm(apply(zzz2, 2, quantile, prob=0.975)),
                              pnorm(apply(zzz3, 2, quantile, prob=0.975)),
                              pnorm(apply(zzz4, 2, quantile, prob=0.975)),
                              pnorm(apply(zzz5, 2, quantile, prob=0.975))),
                     group_ids = rep(c(species_interest[1], species_interest[2], species_interest[3], 
                                   species_interest[4],species_interest[5]), each = N_pred))
                     

indiv_obs <- tibble(counts = indiv$pos / indiv$n,
                    time = indiv$Sampling_Date,
                    group_ids = indiv$Species) 

pool_obs <- tibble(counts = as.numeric(pool$pos > 0),
                    time = pool$Sampling_Date,
                    group_ids = rep('overall', N_pool)) 

rug_data <- tibble(t = c(pool$Sampling_Date, indiv$Sampling_Date), group_ids = c(rep('overall', N_pool), indiv$Species)) 


out_indiv %>% ggplot() +
  geom_ribbon(aes(x=t, ymin=lower, ymax=uper), fill='blue', alpha=0.25) +
  geom_ribbon(aes(x=t, ymin=lower, ymax=uper), fill='blue', alpha=0.25, data = preds_grp1) +
  geom_line(aes(x=t, y=median), linetype = 3, data = preds_grp1) +
  theme_bw() +
  facet_wrap(~factor(group_ids, levels = c(species_interest, 'overall'))) + 
  geom_line( aes(x=t, y=median), linetype=3) +
  geom_point(data=indiv_obs, aes(x=time, y=counts), shape=4) +
  geom_point(data=pool_obs, aes(x=time, y=counts), shape=17, alpha = .2) +
    ylim(0,1) + geom_rug(aes(t), data = rug_data) + ylab("Prevalence") + 
  xlab("Date") + scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  ggtitle("Data integration procedure for Congo Basin bats")
