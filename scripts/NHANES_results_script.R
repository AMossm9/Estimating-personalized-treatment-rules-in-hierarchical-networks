#Script for reading results of NHANES analysis

#Clear environment
rm(list = ls())

#Check if packages are already installed before installing
if(!require("lme4")) install.packages("lme4")

if(!require("ggplot2")) install.packages("ggplot2")

if(!require("latex2exp")) install.packages("latex2exp")

if(!require("MASS")) install.packages("MASS")

if(!require("dplyr")) install.packages("dplyr")

library(lme4)

library(ggplot2)

library(latex2exp)

library(MASS)

library(dplyr)

#Reading in NHANES results

nhanes_btwn_wthn <- readRDS(here::here("results", "nhanes_intf_btwn_and_wthn.rds"))

nhanes_wthn <- readRDS(here::here("results", "nhanes_intf_wthn.rds"))


blips_btwn_wthn <- nhanes_btwn_wthn[[6]]

spillovers_btwn_wthn <- nhanes_btwn_wthn[[7]]

stage_1_blips_btwn_wthn <- c(mean(blips_btwn_wthn[,1], na.rm=TRUE), mean(blips_btwn_wthn[,2], na.rm=TRUE), mean(blips_btwn_wthn[,3], na.rm=TRUE))

stage_2_blips_btwn_wthn <- c(mean(blips_btwn_wthn[,4], na.rm=TRUE), mean(blips_btwn_wthn[,5], na.rm=TRUE), mean(blips_btwn_wthn[,6], na.rm=TRUE))

tailoring_vars <- c("GMI", "GMI*Alter", "GMI*Housing")

blip_results_btwn_wthn <- data.frame(tailoring_vars, stage_1_blips_btwn_wthn, stage_2_blips_btwn_wthn)

names(blip_results_btwn_wthn) <- c("Tailoring Variable", "Stage 1", "Stage 2")

spillover_vals_btwn_wthn <- c(mean(spillovers_btwn_wthn[,1], na.rm=TRUE), mean(spillovers_btwn_wthn[,2], na.rm=TRUE))

stages_btwn_wthn <- c("Stage 1", "Stage 2")

spillover_results_btwn_wthn <- data.frame(stages_btwn_wthn, spillover_vals_btwn_wthn)

names(spillover_results_btwn_wthn) <- c("Stage", "Spillover Estimate")

write.csv(blip_results_btwn_wthn, here::here("results", "nhanes_blips_btwn_wthn.csv"))

write.csv(spillover_results_btwn_wthn, here::here("results", "nhanes_spillover_btwn_wthn.csv"))

#************************************************************************************************

blips_wthn <- nhanes_wthn[[2]]

spillovers_wthn <- nhanes_wthn[[3]]

stage_1_blips_wthn <- c(mean(blips_wthn[,1], na.rm=TRUE), mean(blips_wthn[,2], na.rm=TRUE), mean(blips_wthn[,3], na.rm=TRUE))

stage_2_blips_wthn <- c(mean(blips_wthn[,4], na.rm=TRUE), mean(blips_wthn[,5], na.rm=TRUE), mean(blips_wthn[,6], na.rm=TRUE))

tailoring_vars <- c("GMI", "GMI*Alter", "GMI*Housing")

blip_results_wthn <- data.frame(tailoring_vars, stage_1_blips_wthn, stage_2_blips_wthn)

names(blip_results_wthn) <- c("Tailoring Variable", "Stage 1", "Stage 2")

spillover_vals_wthn <- c(mean(spillovers_wthn[,1], na.rm=TRUE), mean(spillovers_wthn[,2], na.rm=TRUE))

stages_wthn <- c("Stage 1", "Stage 2")

spillover_results_wthn <- data.frame(stages_wthn, spillover_vals_wthn)

names(spillover_results_wthn) <- c("Stage", "Spillover Estimate")

write.csv(blip_results_wthn, here::here("results", "nhanes_blips_wthn.csv"))

write.csv(spillover_results_wthn, here::here("results", "nhanes_spillover_wthn.csv"))


