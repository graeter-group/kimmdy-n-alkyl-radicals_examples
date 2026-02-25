
# n-alkyl radical examples scripts
See [kimmdy-examples](https://github.com/graeter-group/kimmdy-examples) for other examples of the manuscript.

## Overview

* [simulations](..%2simulations) should contain the MD and KIMMDY simulation runs. Note: This example folder only contains basic run files like gro, yaml etc. A full KIMMDY run will have additional output files.
* [output](..%2output) collects the output files from the scripts (like figures)
* [experimental.py](..%2Fexperimental.py) contains the HAT reaction rate equations from various sources.

### n-alkyl molecules
In [templates](simulations/n-alkyl_examples/templates) are the templates from which the n-alkyl simulation files were build using [prepare_files.py](simulations/n-alkyl_examples/prepare_files.py).
1.  [n_alkyl__extract_rates.py](..%2Fn_alkyl__extract_rates.py) collects the simulated n-alkyl radical rates. 
Note: ([simulations/n-alkyl_examples](..%2simulations/n-alkyl_examples) must be populated through kimmdy runs first)
2.  [fig_3c_and_S2.py](fig_3c_and_S2.py) plots figures 3c and S2.

### Hyperparameter runs
[reduce_trajectory.sh](simulations/hyperparameter_grid/reduce_trajectory.sh) was used to subsample from 4 octyl runs of [octyl](simulations/n-alkyl_examples/systems/500K/octyl).
1. [hyperparameters__extract_rates.py](..%2Fhyperparameters__extract_rates.py) collects the 
simulated 1-octyl rates for different run times and number of barrier predictions. 
Note: ([simulations/hyperparameter_grid](..%2simulations/hyperparameter_grid) must be populated through kimmdy runs first)
2. [fig_3d_and_table_S5.py](fig_3d_and_table_S5.py) plots figures 3d and table S5.
