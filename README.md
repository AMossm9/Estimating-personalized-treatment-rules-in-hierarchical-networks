Reproducibility Scripts for: "Estimating personalized treatment rules in hierarchical networks with
interference using dynamic weighted ordinary least squares"

This repository contains all coding scripts to reproduce the simulation results, figures, and table entries in the manuscript "personalized treatment rules in hierarchical networks with
interference using dynamic weighted ordinary least squares" for the International Journal of Biostatistics.

Analyses were conducted using R version 4.5.2 and package dependencies are managed using the 'renv' package for reproducibility.

Required packages:
- lme4 
- ggplot2 
- latex2exp
- MASS
- dplyr
- nhanesA
- DTRreg
- here

All simulations use fixed random seeds for reproducibility.

Estimating-personalized-treatment-rules-in-hierarchical-networks
│
scripts/
├── sims_one_stage_trt_and_trt_fr_correct.R This script provides working functions for the one-stage analysis in Study 1a with the treatment and treatment-free models correctly specified
├── sims_one_stage_trt_and_trt_fr_correct_main.R This script runs the above script and stores the results
├── sims_one_stage_trt_and_trt_fr_mis.R This script provides working functions for the one-stage analysis in Study 1a with the treatment and treatment-free models mis-specified
├── sims_one_stage_trt_and_trt_fr_mis_main.R This script runs the above script and stores the results
├── sims_one_stage_trt_fr_mis.R This script provides working functions for the one-stage analysis in Study 1a with the treatment-free model mis-specified
├── sims_one_stage_trt_fr_mis_main.R This script runs the above script and stores the results
├── sims_one_stage_trt_mis.R This script provides working functions for the one-stage analysis in Study 1a with the treatment model mis-specified
├── sims_one_stage_trt_mis_main.R This script runs the above script and stores the results
├── sims_study1b_1to5.R This script provides working functions for the one-stage analysis in Study 1b
├── sims_study1b_1to5_main.R This script runs the above script and stores the results
├── sims_study1b_6to10.R This script provides working functions for the one-stage analysis in Study 1b
├── sims_study1b_6to10_main.R This script runs the above script and stores the results
├── sims_study1b_extension_trt_fr_mis.R This script provides working functions for the one-stage analysis in Study 1b in the Supplementary File
├── sims_study1b_extension_trt_fr_mis_main.R This script runs the above script and stores the results
├── sims_study1b_extension_trt_mis.R This script provides working functions for the one-stage analysis in Study 1b in the Supplementary File
├── sims_study1b_extension_trt_mis_main.R This script runs the above script and stores the results
├── sims_two_stage_trt_fr_mis.R This script provides working functions for the two-stage analysis in Study 1a with the treatment-free model mis-specified
├── sims_two_stage_trt_fr_mis_variation.r This script provides working functions for the two-stage analysis in Study 1a with the treatment-free model mis-specified
├── sims_two_stage_trt_mis.R This script provides working functions for the two-stage analysis in Study 1a with the treatment model mis-specified
├── figures_script.R This script produces the figures and results for tables in the manuscript
├── NHANES_interference_analysis.R This script provides the analysis for the NHANES example in Section 4
└── NHANES_results_script.R This script provides results for tables for the NHANES example in Section 4
│
results/
├──
│
renv/
│
.Rprofile
│
.gitignore
|
Estimating-personalized-treatment-rules-in-hierarchical-networks.Rproj
|
renv.lock


How to cite this material:

To use this repository, please cite:
