# FIV_FeLV_Transmission

This repository features code and data used in the manuscript "Transmission of one predicts another: Apathogenic proxies for transmission dynamics of a fatal virus"

The analysis pipeline follows the following steps:
1. Phyloscanner analysis of FIV genomes
2. Exponential random graph model (ERGM) analysis of FIV transmission network
3. FeLV transmission simulation based on FIV model and 3 alternative models
4. Analysis of simulation results


## 1. Phyloscanner analysis of FIV genomes
FIV genomes are available in GenBank. See [Malmberg et al 2019](https://doi.org/10.1098/rspb.2019.1689)

Phyloscanner analysis parameters are provided in "Phyloscanner_networks/Phyloscanner run parameters.csv"


## 2. Exponential random graph modeling
Scripts and associate data for exponential random graph modeling follow the following pipeline steps (scripts can be found in "R Scripts" folder):
1. **Phyloscanner_to_network.R:** converts Phyloscanner output to network object. Uses:
    1. "Phyloscanner_networks/FL FIV_hostRelationshipSummary.csv"
2. **UnwUndERGMs.R:** Exponential random graph modeling of FIV transmission network. Requires the following data:
    1. "Attribute Data/FL FIV covariates.Rdata"
    2. "Attribute Data/FL FIV_panther relatedness.Rdata"
    3. "Attribute Data/FL FIV pairwise dists_logkm.Rdata"
    4. "Attribute Data/FL FIV pairwise overlap_UDOI.Rdata"
    5. "Phyloscanner_networks/FL FIV transmission network.Rdata": Phyloscanner transmission network output from Phyloscanner_to_network.R


## 3. FeLV transmission simulation
FeLV transmission is then simulated on FIV-based networks and three alternative models: random networks, spatial overlap-based (SO) networks, and homogeneous mixing. Main simulation scripts call additional external functions as follows (scripts can be found in "R Scripts" folder):
1. **LHS_design.R:** Latin hypercube sampling to generate simulation parameter sets
    1. "LHS parameter sets.Rdata": Actual parameter sets which can be used to reproduce simulations (required by all simulation scripts)
2. **FIV_full_sims.R:** main simulation script for FIV-based model. Calls external functions:
    1. ***simulate_ergm.R:*** simulates network based on FIV ERGM and population characteristics. Calls data:
        1. "Attribute Data/FeLV Period telemetry_deid.Rdata": de-identified telemetry data from the period of the empirical FeLV outbreak.
        2. "Phyloscanner_networks/FL FIV_bestERGM.Rdata": best FIV ERGM from exponential random graph modeling.
    2. ***trans_sim.R:*** simulations FeLV transmission on network
    3. ***post_process_outbreak_data.R:*** processes simulation results to account for re-spawning process
    4. ***props affected_births included.R:*** processes simulation data to provide proportions of population in each disease category for plotting purposes
    5. ***extract_results.R:*** extracts simulation results of interest
3. **Rand_full_sims.R:** main simulation script for random network model. Calls functions b-d from FIV above, as well as:
    1. ***simulate_randnet.R:*** simulates random network
4. **SO_full_sims.R:** main simulation script for overlap-based network model. Calls functions b-d from FIV above, as well as:
    1. ***simulate_SOnet.R:*** simulates overlap based network. Requires:
        1. "SO network sims_nbinom params.Rdata": negative binomial distribution parameters describing average degree distribution for panther spatial overlap networks during the FeLV period.
5. **FeLV_homogmix.R:** main simulation script for homogeneous mixing model (Gillespie algorithm)


## 4. Analysis of simulation results
Simulation results across all model types are then analyzed using the following (scripts can be found in "R Scripts" folder):
1. **Analyze_simulation_results.R:** main analysis script. Identifies "feasible" parameter sets, performs random forest for parameter importance, and performs GLMM for model type performance. Requires LHS parameter data.
2. **SatScan_sims.R:** preps simulation data for SaTScan analysis, then extracts, processes, and plots results. Note: requires external running of [SaTScan](https://www.satscan.org/).
3. **CuzEd_sims.R:** runs Cuzick-Edwards test for spatial clustering with simulation data. 
4. **RF_cbeta_figure.R:** runs random forest analysis to generate supplementary Figure S10.
