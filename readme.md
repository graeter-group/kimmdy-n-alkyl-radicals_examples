
# n-alkyl radical examples scripts
See [kimmdy-examples](https://github.com/graeter-group/kimmdy-examples) for other examples of the manuscript.

## Overview

* [simulations](..%2simulations) should contain the MD and KIMMDY simulation runs
* [output](..%2output) collects the output files from the scripts (like figures)
* [experimental.py](..%2Fexperimental.py) contains the HAT reaction rate equations from various sources.

### n-alkyl molecules
1. Run [n_alkyl__extract_rates.py](..%2Fn_alkyl__extract_rates.py) to collect the 
simulated n-alkyl radical rates. 
Note: ([simulations/n-alkyl_examples](..%2simulations/n-alkyl_examples) must be populated first)
2. Then use [fig_2c_and_A1a.py](..%2Ffig_2c_and_A1a.py) to plot figures 2c and A1a

### Hyperparameter runs
1. Run [hyperparameters__extract_rates.py](..%2Fhyperparameters__extract_rates.py) to collect the 
simulated 1-octyl rates for different run times and number of barrier predictions. 
Note: ([simulations/hyperparameter_grid](..%2simulations/hyperparameter_grid) must be populated first)
2. Then use [fig_2d_and_A1b.py](..%2Ffig_2d_and_A1b.py) to plot figures 2d and A1b
