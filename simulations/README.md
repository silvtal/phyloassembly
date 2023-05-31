# README.md

Neutral simulations can be replicated with the `dilgrowth_goldford.R` script. The code has since then been updated to include more options and make it more generalizable, faster C++ code for the speed-limiting growth step and easier input format. You can find it as an R package [here](https://github.com/silvtal/dilgrowth).

## About the mechanistic simulations
We implemented a computational model that simulates drift in communities undergoing dilution-growth cycles. At each cycle, the community is first randomly subsampled to mimic experimental dilution. Then, 
the community grows by adding one new individual at each step until the pre-dilution community size is reached. At each growth step, the probability for each population to grow equals its actual relative 
abundance. The process thus models neutral growth along a passage experiment where all individuals have equal fitness.

We also introduced a second modelling scenario ("constrained neutrality") where the initial community populations are initially split into different groups, which then are subjected to the above dilution-growth cycles independently with different pre-set per group community sizes. The output contains each final community or all trajectory communities depending on the settings. The number of transfers and dilution rates employed (12 and 0.008, respectively) were those employed by Goldford _et al._ in their experiments. The initial simulated communities reflected the starting community compositions of their experiments and the community size matched its respective sequencing depths. Thus, this scenario limits growth using the observed average relative abundance of each PCGs in the glucose experiment. This abundance is 71.5% for Node35562 (related to Enterobacteriaceae), 21% for Node27828 (related to Pseudomonadaceae) and 7.5% for a group containing all other OTUs.


## Datasets
Real data : the initial ("original") abundances tables and their metadata maps are in `simulation_parameters` alongside the abundance table for endpoint communities (`table_tr12.txt`).

## Main analysis scripts
# TODO meter cada input mencionado en /simulation_results
Simulated data : in my_neutral_model_glc_CLUS_<date>. There are multiple versions, the most recent one is the good one.

1. my_neutral_model(_parallel).R
   - Main script. Runs the simulation described at the upper section. **#TODO: make it more computationally efficient. HMM**.
   - my_neutral_model_parallel_16.R is a faster version, better parallelised.

2. analysis_launcher.R
   - my_neutral_model_analysis_distr.R --> Two analyses: 1) all the OTUs are considered equal, without taking into consideration their PCG. 2) The three PCGs are analysed separately. Computes faith's PD, shannon's div, Pielou's evenness and the observed richness. Plots ECDF (extra options: Kolmigoroff-Smirnov + plot_boxplots inside the ECDF plot)
   - my_neutral_model_analysis_distr_PERC.R --> Reads all the separate simulations (3 simul types, one per PCG, applying fixed percentages), puts them together and obtains metrics and ECDF. 
   - my_neutral_model_analysis_distr_wo_others.R --> Same as _PERC, but excluding others. Not really useful except for double-checking.

3. boxplots_by_metric.R (alternative version: boxplots.R)
   - graphical comparison of real vs. neutral simul vs. neutral simul w/ constraints (applying fixed percentages). For each metric

4. comparisons_by_metric.R (alt ver comparisons_by_X.R)
   - DTS real vs. simul vs. simul w/ constraints.

5. 2021-11-1_RA_survival/Scripts
   - outlier_id.R prints a list of outliers --> OTUs that grow more or less than expected (= according to the simulations)
   - outlier_id_t_test.R same with T-test instead of permutations. Can't handle stdev==0.
   - RAI_survival_simul.R generates "initial relative abundance" vs. "survival rate" plots. Lineal regression. **#TODO** The "limit" option does not work. The idea is that you can actually see the outliers, but without statistical backing.
   - RAI_survival_real.R same but for real data

6. resultados: PCA_PCoA_<fecha>; scripts en PCA_PCoA_scripts
   - (6 scripts in total, the useful ones are:)
   - PCA_SAMPLES_real_vs_simul.R     (input is an abundance table)
   - PCA_SAMPLES_real_vs_simul_PCG.R (same for PCGs)
   - (los demás sí usan métricas, y hay una flecha para cada métrica)
   - (Actually not uploaded here "PCA_PCoA_wu_SAMPLES_ONE_TABLE.R" (temp name) was used for the 2022 paper)

7. good_graphs_cordero
   - (correlation graphs using Goyal et al.'s method)
   - (**!!** these ones were just a test using my simulations instead of the gamma distribution data (as seen in the 2022 paper); those are in a different repository)
   - real_transfer_explorer.R quick check to see for which replicates we have data for at least 5 time points
   - n1_my_neutral_model_parallel.R same as (1.)
   - n2_correlations.R get Pearson correlations for each simulated OTU pair, and each real OTU pair. I use all the available time points available at 8-12, without considering the replicate
   - n3_networks.R picks those correlations, applies a filter to check which correlations are actually relevant. We check with a p-value obtained with a permutations approach (similar to SparCC, Goyal et al.). ORIGINALLY, the "permutations" are actually correlations computed from random values picked from a gamma distribution of abundances of each OTU. IN THIS CASE however, for each OTU we take 800 abundances simulated with (1.) and obtain the Pearson correlations, then we use the distribution of those 800 simuls to obtain a p value

8. plot_lindo.R
   - cute plot explaining how the simulations work


### Functions for the scripts above
   - functions_for_neutral_modelling.R 
   - my_functions.R 

### Additional scripts

   - fix_reports.R y fix_reports_simul.R : generate "reports" for the simulations, which include the rank (abundance) of all OTUs + their survival rate + their initial relative abundance
   - my_neutral_model_perc_applier.R     : puts together simulations for separate PCGs
   - final_samples_PCG_perc_explorer.R   : generates data/map_final_perc.csv
