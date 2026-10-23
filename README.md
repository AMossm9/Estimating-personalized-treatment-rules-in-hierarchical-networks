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
├── sims_study1b_1to5.R This script provides working functions for the one-stage analysis in Study 1b when the treatment model is correctly specified
├── sims_study1b_1to5_main.R This script runs the above script and stores the results
├── sims_study1b_6to10.R This script provides working functions for the one-stage analysis in Study 1b when the treatment model is mis-specified
├── sims_study1b_6to10_main.R This script runs the above script and stores the results
├── sims_study1b_extension_trt_fr_mis.R This script provides working functions for the one-stage analysis in Study 1b when two of the ego's and alters' components in the treatment-free model are mis-specified and the treatment model is correct (Supplementary) 
├── sims_study1b_extension_trt_fr_mis_main.R This script runs the above script and stores the results
├── sims_study1b_extension_trt_mis.R This script provides working functions for the one-stage analysis in Study 1b when two of the ego's and alters' components in the treatment-free model are mis-specified and the treatment model is mis-specified (Supplementary)
├── sims_study1b_extension_trt_mis_main.R This script runs the above script and stores the results
├── sims_two_stage_trt_fr_mis.R This script provides working functions for the two-stage analysis in Study 1a with the treatment-free model mis-specified
├── sims_two_stage_trt_fr_mis_variation.r This script provides working functions for the two-stage analysis in Study 1a with the treatment-free model mis-specified through the cluster-level variable
├── sims_two_stage_trt_mis.R This script provides working functions for the two-stage analysis in Study 1a with the treatment model mis-specified
├── figures_script.R This script produces the figures and results for tables in the manuscript
├── NHANES_interference_analysis.R This script provides the analysis for the NHANES example in Section 4
└── NHANES_results_script.R This script provides results for tables for the NHANES example in Section 4
│
results/
├── nhanes_intf_btwn_and_wthn.rds This script contains between- and within-group interference results from the NHANES_interference_analysis.R script, read in the NHANES_results_script.R
├── nhanes_intf_wthn.rds This script contains within-group only interference results from the NHANES_interference_analysis.R script, read in the NHANES_results_script.R
├── one_stage_trt_and_trt_fr_correct.rds This script contains results from running sims_one_stage_trt_and_trt_fr_correct.R, read in the figures_script.R
├── one_stage_trt_and_trt_fr_mis.rds This script contains results from running sims_one_stage_trt_and_trt_fr_mis.R, read in the figures_script.R
├── one_stage_trt_fr_mis.rds This script contains results from running sims_one_stage_trt_fr_mis.R, read in the figures_script.R
├── one_stage_trt_mis.rds This script contains results from running sims_one_stage_trt_mis.R, read in the figures_script.R
├── sims_two_stage_trt_fr_mis_variation.rds This script contains results from running sims_two_stage_trt_fr_mis_variation.r, read in the figures_script.R
├── sims_two_stage_trt_fr_mis.rds This script contains results from running sims_two_stage_trt_fr_mis.R, read in the figures_script.R
├── sims_two_stage_trt_mis.rds This script contains results from running sims_two_stage_trt_mis.R, read in the figures_script.R
├── study1b_1to5.rds This script contains results from running sims_study1b_1to5.R, read in the figures_script.R
├── study1b_6to10.rds This script contains results from running sims_study1b_6to10.R, read in the figures_script.R
├── study1b_extension_trt_fr_mis.rds This script contains results from running sims_study1b_extension_trt_fr_mis.R, read in the figures_script.R
├── study1b_extension_trt_mis.rds This script contains results from running sims_study1b_extension_trt_mis.R, read in the figures_script.R
└── All other files are the results corresponding to the above .rds files run
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
