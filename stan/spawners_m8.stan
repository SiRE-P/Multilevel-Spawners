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
  // group distributions (years)
  real<lower=0> arrival_mu;         
  real<lower=0> arrival_sigma;           // variation across years
  
  real      arrival_spread_mu;     // stdev entry
  real<lower=0> arrival_spread_sigma;     
  
  real      exit_spread_mu;      // stdev exit
  real<lower=0> exit_spread_sigma;  
  
  real      exit_lag_mu;  // lag entry to exit
  real<lower=0> exit_lag_sigma;  

  // simple
  real<lower=1>         live_phi;   // dispersion parameter on live fish sampling
  
  // instance in group, multi-level
  vector[n_years] log_run;
  vector[n_years] arrival_z;
  vector[n_years] arrival_spread_z;
  vector[n_years] exit_spread_z;
  vector[n_years] exit_lag_z;
}

transformed parameters{
  //non-centered priors
  vector[n_years] arrival;
  vector[n_years] arrival_spread;
  vector[n_years] exit_spread;
  vector[n_years] exit_lag;
  
  
  for (y in 1:n_years){
    arrival[y] = arrival_mu + arrival_sigma * arrival_z[y];
    arrival_spread[y] = exp(arrival_spread_mu + arrival_spread_sigma * arrival_spread_z[y]) + 1;
    exit_spread[y] = exp(exit_spread_mu + exit_spread_sigma * exit_spread_z[y]) + 1;
    exit_lag[y] = exp(exit_lag_mu + exit_lag_sigma * exit_lag_z[y]) + 1;
  }

  
  //process model
  // live
  vector[n_obs] live_mu; // predictions from model
  for(i in 1:n_obs){
    real mean_arrival = arrival[year[i]];
    real entered = normal_cdf(day[i], mean_arrival, arrival_spread[year[i]]);  
    real exited  = normal_cdf(day[i], mean_arrival + exit_lag[year[i]], exit_spread[year[i]]);
    real log_p    = log(fmin(fmax(entered * (1 - exited), 1e-5), 1.0));
    real log_mu   = log_run[year[i]] + log_p;

    live_mu[i]    = exp(log_mu); 
    }
}

model {
  // local variables
  // group prior PDDs
  
  arrival_mu ~ normal(priors[2,1], priors[2,2]); 
  arrival_sigma ~ lognormal(log(priors[3,1]), priors[3,2]);                
  arrival_z ~ normal(0, 1);             
  
  arrival_spread_mu  ~ normal(priors[4,1], priors[4,2]);  // 
  arrival_spread_sigma ~ lognormal(log(priors[5,1]), priors[5,2]);
  arrival_spread_z ~ normal(0, 1);
  
  exit_spread_mu   ~ normal(priors[4,1], priors[4,2]);  // 
  exit_spread_sigma ~ lognormal(log(priors[5,1]), priors[5,2]);
  exit_spread_z ~ normal(0, 1);
  
  exit_lag_mu ~ normal(priors[6,1], priors[6, 2]);
  exit_lag_sigma ~ lognormal(log(priors[5,1]), priors[5,2]);
  exit_lag_z ~ normal(0, 1);
  
  // simple PDD
  log_run ~ normal(priors[1,1], priors[1,2]);  // mean  run log
  log(live_phi)   ~ normal(priors[7,1], priors[7,2]);
  
  
  //likelihood
  target += neg_binomial_2_lpmf(live_counts | live_mu, live_phi);
}

generated quantities{
  vector[n_obs] log_lik;
  for(i in 1:n_obs)
  log_lik[i] = neg_binomial_2_lpmf(live_counts[i] | live_mu[i], live_phi);
}
