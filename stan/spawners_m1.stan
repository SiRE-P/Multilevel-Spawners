data {
  int                    n_priors;        // n priors
  matrix[n_priors,2]     priors;          // rows: priors; cols: p1, p2
  int                    n_years;         // n years
  int                    n_obs;          // number of spawner count obs
  
  int<lower=1, upper=n_years>     year[n_obs];
  int                    day[n_obs];
  int                    live_counts[n_obs]; //live spawner counts
}
parameters {

  // simple
  real<lower=0> arrival;         
  real<lower=0>      arrival_spread;     // stdev entry
  real<lower=0>      exit_spread;      // stdev exit
  real     exit_lag_raw;  // lag entry to exit
  real<lower=0>         live_phi;   // dispersion parameter on live fish sampling
  
  // instance in group, multi-level
  vector[n_years] log_run;
}

transformed parameters{
  real exit_lag;
  exit_lag = exp(exit_lag_raw) + 1;
  
  //process model
  // live
  vector[n_obs] live_mu; // predictions from model
  for(i in 1:n_obs){
    real mean_arrival = arrival;
    real entered = normal_cdf(day[i], mean_arrival, arrival_spread);  
    real exited  = normal_cdf(day[i], mean_arrival + exit_lag, exit_spread);
    real log_p    = log(fmin(fmax(entered * (1 - exited), 1e-5), 1.0));
    real log_mu   = log_run[year[i]] + log_p;

    live_mu[i]    = exp(log_mu);  
    }
}

model {
  // local variables
  // group prior PDDs
  
  arrival ~ normal(priors[2,1], priors[2,2]); 
  
  log(arrival_spread) ~ normal(priors[4,1], priors[4,2]);  // 

  log(exit_spread)  ~ normal(priors[4,1], priors[4,2]);  // 
  
  // simple PDD
  log_run ~ normal(priors[1,1], priors[1,2]);  // mean  run log
  exit_lag_raw  ~ normal(priors[6,1], priors[6,2]);  // same all years
  log(live_phi)   ~ normal(priors[7,1], priors[7,2]);
  
  
  //likelihood
  target += neg_binomial_2_lpmf(live_counts | live_mu, live_phi);
}

generated quantities{
  vector[n_obs] log_lik;
  for(i in 1:n_obs)
  log_lik[i] = neg_binomial_2_lpmf(live_counts[i] | live_mu[i], live_phi);
}

