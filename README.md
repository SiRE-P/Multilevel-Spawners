# Multilevel Spawners  
## Objective  
Analyze surveys that count spawning salmon: fit run (spawner abundance), arrival_timing, arrival_spread, and exit_lag, exit_spread.  
## Case  
Wild Sockeye Salmon spawning in the Okanagan River above Osoyoos Lake (natal lake); 9 to 19 surveys/year, 2000 to 2024, total 291 surveys.  
![map](/figures/Okanagan_to_Osoyoo_Lakes.png)  
## Models  
### Simulations  
Within a year, spawners arrive according to a cumulative normal distribution: annual values for arrival timing and spread. They exit similarly: annual values of exit_lag and exit_spread. Spawners present during quasi-random survey dates are observed with sampling error. A suite of prior probability distributions is used for random draws of parameters for simulation. The simulation results lead to data and priors used to fit the simulated data to the model, then resulting fit of surveys, runs, and annual parameters are compared to simulated.  
Code in GitHub: Multilevel-Spawners/r/  
### Fit model to data  
Non-linear, multi-level Bayesian regression (Stan, rstan package) where parameters vary by year. Multi-level means estimates (posterior distributions) of annual values for a parameters are in a group with a fitted distribution for those values. This improves accuracy of otherwise poorly estimated annual parameters.  
Code in GitHub: Multilevel-Spawners/stan  
## Primary Publication  
Drafts in GitHub: Multilevel-Spawners/quarto/  
### Quarto template for primary  
see [https://github.com/quarto-journals/elsevier](https://github.com/quarto-journals/elsevier)
Installation for an existing Quarto project or document: 
- From the quarto project or document directory, run the following command to install this format:
    - quarto add quarto-journals/elsevier
- In your document yaml:  
format:  
..pdf: default  
..elsevier-pdf:  
....keep-tex: true  
see preview of rendered template at [https://quarto-journals.github.io/elsevier/](https://quarto-journals.github.io/elsevier/).
