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
  // MVN timing parameters
  vector[4] timing_mu;
  cholesky_factor_corr[4] L_Omega;
  vector<lower=0>[4] timing_sigma;
  matrix[n_years, 4] timing_raw;
  
  vector[n_years] log_run;
  real<lower=1>         live_phi;   // dispersion parameter on live fish sampling
  
}

transformed parameters{
  array[n_years] vector[4] timing;
  matrix[4, 4] L_Sigma = diag_pre_multiply(timing_sigma, L_Omega);
  vector[n_years] arrival;
  vector[n_years] arrival_spread;
  vector[n_years] exit_lag;
  vector[n_years] exit_spread;
  
  for (y in 1:n_years){
    vector[4] timing_z = timing_raw[y]';  // convert row to vector
    timing[y] = timing_mu + L_Sigma * timing_z;
    arrival[y] = timing[y,1];
    arrival_spread[y] = exp(timing[y,2]) + 1;
    exit_spread[y] = exp(timing[y,3]) + 1;
    exit_lag[y] = exp(timing[y,4]) + 1;
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
  timing_mu[1] ~ normal(priors[2,1], priors[2,2]);  // arrival_mu
  timing_mu[2] ~ normal(priors[4,1], priors[4,2]);  // arrival_spread_mu
  timing_mu[3] ~ normal(priors[6,1], priors[6,2]);  // exit_lag_mu
  timing_mu[4] ~ normal(priors[4,1], priors[4,2]);  // exit_spread_mu
  
  timing_sigma[1] ~ lognormal(log(priors[3,1]), priors[3,2]);  // arrival_sigma
  timing_sigma[2] ~ lognormal(log(priors[5,1]), priors[5,2]);  // arrival_spread_sigma
  timing_sigma[3] ~ lognormal(log(priors[5,1]), priors[5,2]);  // exit_lag_sigma
  timing_sigma[4] ~ lognormal(log(priors[5,1]), priors[5,2]);  // exit_spread_sigma
  
  L_Omega ~ lkj_corr_cholesky(2);
  
  to_vector(timing_raw) ~ std_normal();
  
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
  
  corr_matrix[4] Omega;
  cov_matrix[4] Sigma;
  Omega = multiply_lower_tri_self_transpose(L_Omega);  // reconstruct correlation matrix
  Sigma = quad_form_diag(Omega, timing_sigma);         // full covariance matrix
}
