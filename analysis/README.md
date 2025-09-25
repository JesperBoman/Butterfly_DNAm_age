**Varying seed clock pipeline**

#Preparing data scripts:

1. Comb_data.sh
2. add_E2_data.sh

See example data structure in: example_input_data.txt

#Varying seed clock
This pipeline randomly splits training+test data based on a particular seed and constructs a clock for each seed. We do this to evaluate the effect of a particular training+test split on predictiveness of the model. In a sense, we then get sets of clocks. 

master.clocksmith.sh runs smither.sh which in turn runs fancier_varying_seed_clock.R

It is fancier_varying_seed_clock.R that is the script that is training the penalized regression model. The master.clocksmith and smither are used to parallelize the process without sending all jobs at once, to avoid sending too many concurrent jobs to a cluster.
