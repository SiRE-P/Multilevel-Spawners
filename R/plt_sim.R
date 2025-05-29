library(tidyverse)
library(rstan)
library(tidybayes)
set.seed(123)

log_run_mu <- 10
arrival_mu <- 35
arrival_spread_mu <- 5
exit_spread_mu <- 5
exit_lag <- 10

years <- 10
N_total <- round(rlnorm(years, log_run_mu, 1))

arrival <- rnorm(years, arrival_mu, 5) 
arrival_spread <- rlnorm(years, log(arrival_spread_mu), 0.3)
exit_spread <- rlnorm(years, log(exit_spread_mu), 0.3)
exit_lag <- rlnorm(years, log(exit_lag), 0.3)

output <- sapply(1:years, function(i){
  x <- 1:100
  entered <- pnorm(x, arrival[i], arrival_spread[i])
  exited <- pnorm(x, arrival[i] + exit_lag[i], exit_spread[i])
  
  present <- entered * (1-exited) * N_total[i]
})

matplot(output, type = 'l')

full_surveyed.df <- data.frame(N = round(c(output)), sampled = rnbinom(n = length(output), mu = output, size = 10), year = rep(1:years, each = 100), day = 1:100)

sampled.df <- full_surveyed.df %>% 
  group_by(year) %>% 
  sample_frac(size = 0.9)

ggplot(sampled.df, aes(x = day, y = sampled))+
  geom_point()+
  facet_wrap(~year)

sampled.df2 <- sampled.df %>% 
  filter(sampled >0)

priors <- data.frame(prior = c("log_runs_mu", "arrival_mu", "arrival_sigma", "spread","spread_sigma","residence", "count_dispersion"),
                     v1 = c(10,40, 0,log(8-1), 0,log(11-1), 1), 
                     v2 = c(1,5,10,0.2, 5, 0.3, 0.2))

sp_dat = list(n_priors = nrow(priors), 
              priors = data.matrix(priors[,-1]),
              n_years = max(sampled.df2$year),
              year = as.numeric(factor(sampled.df2$year)),
              day = sampled.df2$day,
              n_obs = nrow(sampled.df2),
              live_counts = sampled.df2$sampled, 
              live_phi = 10,
              arrival_mu = arrival_mu)


mod <- stan_model("./stan/spawners_m8.stan")
fit <- sampling(mod, data = sp_dat, chains = 4, cores = 4, iter = 500, seed = 7, init = replicate(4, list(arrival_mu = 35), simplify = FALSE))

traceplot(fit)

post <- extract(fit)

spread_draws(fit, log_run[year]) %>% 
  filter(.chain %in% c(3:4)) %>% 
  ggplot(aes(x = year, y = exp(log_run)))+
  stat_pointinterval()+
  geom_point(data = data.frame(year = 1:years, log_run = log(N_total)), color = 2)

spread_draws(fit, arrival[year]) %>% 
  filter(.chain %in% c(3:4)) %>% 
  ggplot(aes(x = year, y = arrival))+
  stat_pointinterval()+
  geom_point(data = data.frame(year = 1:years, arrival = arrival), color = 2)

spread_draws(fit, arrival_spread[year]) %>% 
  filter(.chain %in% c(3:4)) %>% 
  ggplot(aes(x = year, y = arrival_spread))+
  stat_pointinterval()+
  geom_point(data = data.frame(year = 1:years, arrival_spread = arrival_spread), color = 2)

spread_draws(fit, exit_lag[year]) %>% 
  filter(.chain %in% c(3:4)) %>% 
  ggplot(aes(x = year, y = exit_lag))+
  stat_pointinterval()+
  geom_point(data = data.frame(year = 1:years, exit_lag = exit_lag), color = 2)

spread_draws(fit, exit_spread[year]) %>% 
  filter(.chain %in% c(3:4)) %>% 
  ggplot(aes(x = year, y = exit_spread))+
  stat_pointinterval()+
  geom_point(data = data.frame(year = 1:years, exit_spread = exit_spread), color = 2)


spawn_curves.df <- data.frame()
for(j in 1:ncol(post$log_run)){
  live <- sapply(1:100, FUN = function(x) {
    entered <- pnorm(x, post$arrival[,j], post$arrival_spread[,j])
    exited <- pnorm(x, post$arrival[,j] + post$exit_lag[j], post$exit_spread[,j])
    
    exp(post$log_run[,j] + log(entered - (entered * exited)))
  })
  
  spawn_curves.df <- bind_rows(spawn_curves.df, data.frame(year = j, 
                                                           day = 1:100,
                                                           fish = apply(live, 2, median),
                                                           l95 = apply(live, 2, quantile, probs = 0.025),
                                                           u95 = apply(live, 2, quantile, probs = 0.975)))
}

ggplot(spawn_curves.df, aes(x = day, y = fish))+
  geom_line()+
  geom_ribbon(aes(ymin = l95, ymax = u95), color = NA, alpha = 0.2)+
  facet_wrap(~year, scales = "free_y")+
  geom_point(data = sampled.df2, aes(x = day, y = sampled), size = 0.8)+
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k"), name = "Spawners present (thousands)", expand = expansion(mult = c(0, 0.02)))+
  theme_bw()+
  theme(strip.background = element_rect(color="NA", fill="NA"))+
  scale_x_continuous(labels = ~ format(as.Date(.x, origin = "2023-12-31"), "%b"), 
                     breaks = as.numeric(as.Date(c("2024-09-01", "2024-10-01", "2024-11-01")) - as.Date("2023-12-31")), name = "")+
  theme(panel.grid.minor = element_blank())
