#Script for Generating Figures of One- and Two-Stage Simulations

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

#One Stage: Trt Mis-specified

one_stage_trt_mis <- readRDS(here::here("results", "one_stage_trt_mis.rds"))

names(one_stage_trt_mis) <- c("AOpts_one_trt_mis", "YOpts_one_trt_mis", "Regret_one_trt_mis", "Psi_M1_one_trt_mis", "Psi_M2_one_trt_mis", "Psi_M3_one_trt_mis", "Xsi_M1_one_trt_mis", "Xsi_M2_one_trt_mis", "Xsi_M3_one_trt_mis", "Min_alters_one_trt_mis", "Max_alters_one_trt_mis", "Mean_alters_one_trt_mis", "Median_alters_one_trt_mis")

AOpts.df <- one_stage_trt_mis$AOpts_one_trt_mis
YOpts.df <- one_stage_trt_mis$YOpts_one_trt_mis
Regret.df <- one_stage_trt_mis$Regret_one_trt_mis
Psi.M1 <- one_stage_trt_mis$Psi_M1_one_trt_mis
Psi.M2 <- one_stage_trt_mis$Psi_M2_one_trt_mis
Psi.M3 <- one_stage_trt_mis$Psi_M3_one_trt_mis
Xsi.M1 <- one_stage_trt_mis$Xsi_M1_one_trt_mis
Xsi.M2 <- one_stage_trt_mis$Xsi_M2_one_trt_mis
Xsi.M3 <- one_stage_trt_mis$Xsi_M3_one_trt_mis

#Provides numerical, tabular summary of optimal treatment rates and regrets across methods
plot_trt_rates <- na.omit(AOpts.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(TrtRate),
         mean = mean(TrtRate),
         sd = sd(TrtRate),
         empse = sqrt((sum((TrtRate-mean(TrtRate))^2))/(n-1)),
         .groups="keep")

plot_trt_rates

plot_outcomes <- na.omit(YOpts.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(Outcome),
         mean = mean(Outcome),
         sd = sd(Outcome),
         empse = sqrt((sum((Outcome-mean(Outcome))^2))/(n-1)),
         .groups="keep")

plot_outcomes

plot_regrets <- na.omit(Regret.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(Regret),
         mean = mean(Regret),
         sd = sd(Regret),
         empse = sqrt((sum((Regret-mean(Regret))^2))/(n-1)),
         .groups="keep")

plot_regrets

write.csv(plot_trt_rates, here::here("results", "sims_one_trt_mis_trtRate_tbl.csv"))

write.csv(plot_outcomes, here::here("results", "sims_one_trt_mis_outcome_tbl.csv"))

write.csv(plot_regrets, here::here("results", "sims_one_trt_mis_regret_tbl.csv"))

AOpts.df$phys.grps <- as.factor(AOpts.df$phys.grps)
YOpts.df$phys.grps <- as.factor(YOpts.df$phys.grps)

Psi.M1$phys.grps <- as.factor(Psi.M1$phys.grps)
Psi.M1$Method <- rep("ME2", nrow(Psi.M1))

Psi.M1$phys.grps <- as.factor(Psi.M1$phys.grps)
Psi.M1$Method <- rep("ME2", nrow(Psi.M1))

Psi.M2$phys.grps <- as.factor(Psi.M2$phys.grps)
Psi.M2$Method <- rep("ME3", nrow(Psi.M2))

Psi.M3$phys.grps <- as.factor(Psi.M3$phys.grps)
Psi.M3$Method <- rep("ME4", nrow(Psi.M3))

Xsi.M1$phys.grps <- as.factor(Xsi.M1$phys.grps)
Xsi.M1$Method <- rep("ME2", nrow(Xsi.M1))

Xsi.M2$phys.grps <- as.factor(Xsi.M2$phys.grps)
Xsi.M2$Method <- rep("ME3", nrow(Xsi.M2))

Xsi.M3$phys.grps <- as.factor(Xsi.M3$phys.grps)
Xsi.M3$Method <- rep("ME4", nrow(Xsi.M3))


Psi.0.df <- rbind(Psi.M1[Psi.M1$PsiParameter=="psi.0",], Psi.M2[Psi.M2$PsiParameter=="psi.0",], Psi.M3[Psi.M3$PsiParameter=="psi.0",])
Psi.1.df <- rbind(Psi.M1[Psi.M1$PsiParameter=="psi.1",], Psi.M2[Psi.M2$PsiParameter=="psi.1",], Psi.M3[Psi.M3$PsiParameter=="psi.1",])
Psi.2.df <- rbind(Psi.M1[Psi.M1$PsiParameter=="psi.2",], Psi.M2[Psi.M2$PsiParameter=="psi.2",], Psi.M3[Psi.M3$PsiParameter=="psi.2",])
Psi.3.df <- rbind(Psi.M1[Psi.M1$PsiParameter=="psi.3",], Psi.M2[Psi.M2$PsiParameter=="psi.3",], Psi.M3[Psi.M3$PsiParameter=="psi.3",])

Xsi.0.df <- rbind(Xsi.M1[Xsi.M1$XsiParameter=="xsi.0",], Xsi.M2[Xsi.M2$XsiParameter=="xsi.0",], Xsi.M3[Xsi.M3$XsiParameter=="xsi.0",])
Xsi.1.df <- rbind(Xsi.M1[Xsi.M1$XsiParameter=="xsi.1",], Xsi.M2[Xsi.M2$XsiParameter=="xsi.1",], Xsi.M3[Xsi.M3$XsiParameter=="xsi.1",])

bp.trt <- ggplot(AOpts.df, aes(x=phys.grps, y=TrtRate, fill=Method)) + geom_boxplot() + labs(y="Optimal treatment rate", x="J") + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(text = element_text(size = 14))
bp.out <- ggplot(YOpts.df, aes(x=phys.grps, y=Outcome, fill=Method)) + geom_boxplot() + labs(y="Optimal outcome", x="J") + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(text = element_text(size = 14))
bp.psi0 <- ggplot(Psi.0.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\psi_{0}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi1 <- ggplot(Psi.1.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\psi_{1}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi2 <- ggplot(Psi.2.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\psi_{2}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=-1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi3 <- ggplot(Psi.3.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\psi_{3}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=-1.5, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.xsi0 <- ggplot(Xsi.0.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\xi_{0}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.xsi1 <- ggplot(Xsi.1.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\xi_{1}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))

ggsave(here::here("results", "sims_one_trt_mis_trtRate_bp.pdf"), plot = bp.trt)
ggsave(here::here("results", "sims_one_trt_mis_outcome_bp.pdf"), plot = bp.out)
ggsave(here::here("results", "sims_one_trt_mis_psi0_bp.pdf"), plot = bp.psi0)
ggsave(here::here("results", "sims_one_trt_mis_psi1_bp.pdf"), plot = bp.psi1)
ggsave(here::here("results", "sims_one_trt_mis_psi2_bp.pdf"), plot = bp.psi2)
ggsave(here::here("results", "sims_one_trt_mis_psi3_bp.pdf"), plot = bp.psi3)
ggsave(here::here("results", "sims_one_trt_mis_xsi0_bp.pdf"), plot = bp.xsi0)
ggsave(here::here("results", "sims_one_trt_mis_xsi1_bp.pdf"), plot = bp.xsi1)

ggsave(here::here("results", "sims_one_trt_mis_trtRate_bp.jpg"), plot = bp.trt)
ggsave(here::here("results", "sims_one_trt_mis_outcome_bp.jpg"), plot = bp.out)
ggsave(here::here("results", "sims_one_trt_mis_psi0_bp.jpg"), plot = bp.psi0)
ggsave(here::here("results", "sims_one_trt_mis_psi1_bp.jpg"), plot = bp.psi1)
ggsave(here::here("results", "sims_one_trt_mis_psi2_bp.jpg"), plot = bp.psi2)
ggsave(here::here("results", "sims_one_trt_mis_psi3_bp.jpg"), plot = bp.psi3)
ggsave(here::here("results", "sims_one_trt_mis_xsi0_bp.jpg"), plot = bp.xsi0)
ggsave(here::here("results", "sims_one_trt_mis_xsi1_bp.jpg"), plot = bp.xsi1)

#One Stage: Trt Free Mis-specified

one_stage_trt_fr_mis <- readRDS(here::here("results", "one_stage_trt_fr_mis.rds"))

names(one_stage_trt_fr_mis) <- c("AOpts_one_trt_fr_mis", "YOpts_one_trt_fr_mis", "Regret_one_trt_fr_mis", "Psi_M1_one_trt_fr_mis", "Psi_M2_one_trt_fr_mis", "Psi_M3_one_trt_fr_mis", "Xsi_M1_one_trt_fr_mis", "Xsi_M2_one_trt_fr_mis", "Xsi_M3_one_trt_fr_mis", "Min_alters_one_trt_fr_mis", "Max_alters_one_trt_fr_mis", "Mean_alters_one_trt_fr_mis", "Median_alters_one_trt_fr_mis")

AOpts.df <- one_stage_trt_fr_mis$AOpts_one_trt_fr_mis
YOpts.df <- one_stage_trt_fr_mis$YOpts_one_trt_fr_mis
Regret.df <- one_stage_trt_fr_mis$Regret_one_trt_fr_mis
Psi.M1 <- one_stage_trt_fr_mis$Psi_M1_one_trt_fr_mis
Psi.M2 <- one_stage_trt_fr_mis$Psi_M2_one_trt_fr_mis
Psi.M3 <- one_stage_trt_fr_mis$Psi_M3_one_trt_fr_mis
Xsi.M1 <- one_stage_trt_fr_mis$Xsi_M1_one_trt_fr_mis
Xsi.M2 <- one_stage_trt_fr_mis$Xsi_M2_one_trt_fr_mis
Xsi.M3 <- one_stage_trt_fr_mis$Xsi_M3_one_trt_fr_mis

#Provides numerical, tabular summary of optimal treatment rates and regrets across methods
plot_trt_rates <- na.omit(AOpts.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(TrtRate),
         mean = mean(TrtRate),
         sd = sd(TrtRate),
         empse = sqrt((sum((TrtRate-mean(TrtRate))^2))/(n-1)),
         .groups="keep")

plot_trt_rates

plot_outcomes <- na.omit(YOpts.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(Outcome),
         mean = mean(Outcome),
         sd = sd(Outcome),
         empse = sqrt((sum((Outcome-mean(Outcome))^2))/(n-1)),
         .groups="keep")

plot_outcomes

plot_regrets <- na.omit(Regret.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(Regret),
         mean = mean(Regret),
         sd = sd(Regret),
         empse = sqrt((sum((Regret-mean(Regret))^2))/(n-1)),
         .groups="keep")

plot_regrets

write.csv(plot_trt_rates, here::here("results", "sims_one_trt_fr_mis_trtRate_tbl.csv"))

write.csv(plot_outcomes, here::here("results", "sims_one_trt_fr_mis_outcome_tbl.csv"))

write.csv(plot_regrets, here::here("results", "sims_one_trt_fr_mis_regret_tbl.csv"))

AOpts.df$phys.grps <- as.factor(AOpts.df$phys.grps)
YOpts.df$phys.grps <- as.factor(YOpts.df$phys.grps)

Psi.M1$phys.grps <- as.factor(Psi.M1$phys.grps)
Psi.M1$Method <- rep("ME2", nrow(Psi.M1))

Psi.M2$phys.grps <- as.factor(Psi.M2$phys.grps)
Psi.M2$Method <- rep("ME3", nrow(Psi.M2))

Psi.M3$phys.grps <- as.factor(Psi.M3$phys.grps)
Psi.M3$Method <- rep("ME4", nrow(Psi.M3))

Xsi.M1$phys.grps <- as.factor(Xsi.M1$phys.grps)
Xsi.M1$Method <- rep("ME2", nrow(Xsi.M1))

Xsi.M2$phys.grps <- as.factor(Xsi.M2$phys.grps)
Xsi.M2$Method <- rep("ME3", nrow(Xsi.M2))

Xsi.M3$phys.grps <- as.factor(Xsi.M3$phys.grps)
Xsi.M3$Method <- rep("ME4", nrow(Xsi.M3))


Psi.0.df <- rbind(Psi.M1[Psi.M1$PsiParameter=="psi.0",], Psi.M2[Psi.M2$PsiParameter=="psi.0",], Psi.M3[Psi.M3$PsiParameter=="psi.0",])
Psi.1.df <- rbind(Psi.M1[Psi.M1$PsiParameter=="psi.1",], Psi.M2[Psi.M2$PsiParameter=="psi.1",], Psi.M3[Psi.M3$PsiParameter=="psi.1",])
Psi.2.df <- rbind(Psi.M1[Psi.M1$PsiParameter=="psi.2",], Psi.M2[Psi.M2$PsiParameter=="psi.2",], Psi.M3[Psi.M3$PsiParameter=="psi.2",])
Psi.3.df <- rbind(Psi.M1[Psi.M1$PsiParameter=="psi.3",], Psi.M2[Psi.M2$PsiParameter=="psi.3",], Psi.M3[Psi.M3$PsiParameter=="psi.3",])

Xsi.0.df <- rbind(Xsi.M1[Xsi.M1$XsiParameter=="xsi.0",], Xsi.M2[Xsi.M2$XsiParameter=="xsi.0",], Xsi.M3[Xsi.M3$XsiParameter=="xsi.0",])
Xsi.1.df <- rbind(Xsi.M1[Xsi.M1$XsiParameter=="xsi.1",], Xsi.M2[Xsi.M2$XsiParameter=="xsi.1",], Xsi.M3[Xsi.M3$XsiParameter=="xsi.1",])

bp.trt <- ggplot(AOpts.df, aes(x=phys.grps, y=TrtRate, fill=Method)) + geom_boxplot() + labs(y="Optimal treatment rate", x="J") + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(text = element_text(size = 14))
bp.out <- ggplot(YOpts.df, aes(x=phys.grps, y=Outcome, fill=Method)) + geom_boxplot() + labs(y="Optimal outcome", x="J") + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(text = element_text(size = 14))
bp.psi0 <- ggplot(Psi.0.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\psi_{0}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi1 <- ggplot(Psi.1.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\psi_{1}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi2 <- ggplot(Psi.2.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\psi_{2}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=-1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi3 <- ggplot(Psi.3.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\psi_{3}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=-1.5, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.xsi0 <- ggplot(Xsi.0.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\xi_{0}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.xsi1 <- ggplot(Xsi.1.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\xi_{1}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))

ggsave(here::here("results", "sims_one_trt_fr_mis_trtRate_bp.pdf"), plot = bp.trt)
ggsave(here::here("results", "sims_one_trt_fr_mis_outcome_bp.pdf"), plot = bp.out)
ggsave(here::here("results", "sims_one_trt_fr_mis_psi0_bp.pdf"), plot = bp.psi0)
ggsave(here::here("results", "sims_one_trt_fr_mis_psi1_bp.pdf"), plot = bp.psi1)
ggsave(here::here("results", "sims_one_trt_fr_mis_psi2_bp.pdf"), plot = bp.psi2)
ggsave(here::here("results", "sims_one_trt_fr_mis_psi3_bp.pdf"), plot = bp.psi3)
ggsave(here::here("results", "sims_one_trt_fr_mis_xsi0_bp.pdf"), plot = bp.xsi0)
ggsave(here::here("results", "sims_one_trt_fr_mis_xsi1_bp.pdf"), plot = bp.xsi1)

ggsave(here::here("results", "sims_one_trt_fr_mis_trtRate_bp.jpg"), plot = bp.trt)
ggsave(here::here("results", "sims_one_trt_fr_mis_outcome_bp.jpg"), plot = bp.out)
ggsave(here::here("results", "sims_one_trt_fr_mis_psi0_bp.jpg"), plot = bp.psi0)
ggsave(here::here("results", "sims_one_trt_fr_mis_psi1_bp.jpg"), plot = bp.psi1)
ggsave(here::here("results", "sims_one_trt_fr_mis_psi2_bp.jpg"), plot = bp.psi2)
ggsave(here::here("results", "sims_one_trt_fr_mis_psi3_bp.jpg"), plot = bp.psi3)
ggsave(here::here("results", "sims_one_trt_fr_mis_xsi0_bp.jpg"), plot = bp.xsi0)
ggsave(here::here("results", "sims_one_trt_fr_mis_xsi1_bp.jpg"), plot = bp.xsi1)

#One Stage: Both Trt and Trt Free Models Mis-specified

one_stage_trt_and_trt_fr_mis <- readRDS(here::here("results", "one_stage_trt_and_trt_fr_mis.rds"))

names(one_stage_trt_and_trt_fr_mis) <- c("AOpts_one_trt_and_trt_fr_mis", "YOpts_one_trt_and_trt_fr_mis", "Regret_one_trt_and_trt_fr_mis", "Psi_M1_one_trt_and_trt_fr_mis", "Psi_M2_one_trt_and_trt_fr_mis", "Psi_M3_one_trt_and_trt_fr_mis", "Xsi_M1_one_trt_and_trt_fr_mis", "Xsi_M2_one_trt_and_trt_fr_mis", "Xsi_M3_one_trt_and_trt_fr_mis", "Min_alters_one_trt_and_trt_fr_mis", "Max_alters_one_trt_and_trt_fr_mis", "Mean_alters_one_trt_and_trt_fr_mis", "Median_alters_one_trt_and_trt_fr_mis")

AOpts.df <- one_stage_trt_and_trt_fr_mis$AOpts_one_trt_and_trt_fr_mis
YOpts.df <- one_stage_trt_and_trt_fr_mis$YOpts_one_trt_and_trt_fr_mis
Regret.df <- one_stage_trt_and_trt_fr_mis$Regret_one_trt_and_trt_fr_mis
Psi.M1 <- one_stage_trt_and_trt_fr_mis$Psi_M1_one_trt_and_trt_fr_mis
Psi.M2 <- one_stage_trt_and_trt_fr_mis$Psi_M2_one_trt_and_trt_fr_mis
Psi.M3 <- one_stage_trt_and_trt_fr_mis$Psi_M3_one_trt_and_trt_fr_mis
Xsi.M1 <- one_stage_trt_and_trt_fr_mis$Xsi_M1_one_trt_and_trt_fr_mis
Xsi.M2 <- one_stage_trt_and_trt_fr_mis$Xsi_M2_one_trt_and_trt_fr_mis
Xsi.M3 <- one_stage_trt_and_trt_fr_mis$Xsi_M3_one_trt_and_trt_fr_mis

#Provides numerical, tabular summary of optimal treatment rates and regrets across methods
plot_trt_rates <- na.omit(AOpts.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(TrtRate),
         mean = mean(TrtRate),
         sd = sd(TrtRate),
         empse = sqrt((sum((TrtRate-mean(TrtRate))^2))/(n-1)),
         .groups="keep")

plot_trt_rates

plot_outcomes <- na.omit(YOpts.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(Outcome),
         mean = mean(Outcome),
         sd = sd(Outcome),
         empse = sqrt((sum((Outcome-mean(Outcome))^2))/(n-1)),
         .groups="keep")

plot_outcomes

plot_regrets <- na.omit(Regret.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(Regret),
         mean = mean(Regret),
         sd = sd(Regret),
         empse = sqrt((sum((Regret-mean(Regret))^2))/(n-1)),
         .groups="keep")

plot_regrets

write.csv(plot_trt_rates, here::here("results", "sims_one_trt_and_trt_fr_mis_trtRate_tbl.csv"))

write.csv(plot_outcomes, here::here("results", "sims_one_trt_and_trt_fr_mis_outcome_tbl.csv"))

write.csv(plot_regrets, here::here("results", "sims_one_trt_and_trt_fr_mis_regret_tbl.csv"))

AOpts.df$phys.grps <- as.factor(AOpts.df$phys.grps)
YOpts.df$phys.grps <- as.factor(YOpts.df$phys.grps)

Psi.M1$phys.grps <- as.factor(Psi.M1$phys.grps)
Psi.M1$Method <- rep("ME2", nrow(Psi.M1))

Psi.M2$phys.grps <- as.factor(Psi.M2$phys.grps)
Psi.M2$Method <- rep("ME3", nrow(Psi.M2))

Psi.M3$phys.grps <- as.factor(Psi.M3$phys.grps)
Psi.M3$Method <- rep("ME4", nrow(Psi.M3))

Xsi.M1$phys.grps <- as.factor(Xsi.M1$phys.grps)
Xsi.M1$Method <- rep("ME2", nrow(Xsi.M1))

Xsi.M2$phys.grps <- as.factor(Xsi.M2$phys.grps)
Xsi.M2$Method <- rep("ME3", nrow(Xsi.M2))

Xsi.M3$phys.grps <- as.factor(Xsi.M3$phys.grps)
Xsi.M3$Method <- rep("ME4", nrow(Xsi.M3))


Psi.0.df <- rbind(Psi.M1[Psi.M1$PsiParameter=="psi.0",], Psi.M2[Psi.M2$PsiParameter=="psi.0",], Psi.M3[Psi.M3$PsiParameter=="psi.0",])
Psi.1.df <- rbind(Psi.M1[Psi.M1$PsiParameter=="psi.1",], Psi.M2[Psi.M2$PsiParameter=="psi.1",], Psi.M3[Psi.M3$PsiParameter=="psi.1",])
Psi.2.df <- rbind(Psi.M1[Psi.M1$PsiParameter=="psi.2",], Psi.M2[Psi.M2$PsiParameter=="psi.2",], Psi.M3[Psi.M3$PsiParameter=="psi.2",])
Psi.3.df <- rbind(Psi.M1[Psi.M1$PsiParameter=="psi.3",], Psi.M2[Psi.M2$PsiParameter=="psi.3",], Psi.M3[Psi.M3$PsiParameter=="psi.3",])

Xsi.0.df <- rbind(Xsi.M1[Xsi.M1$XsiParameter=="xsi.0",], Xsi.M2[Xsi.M2$XsiParameter=="xsi.0",], Xsi.M3[Xsi.M3$XsiParameter=="xsi.0",])
Xsi.1.df <- rbind(Xsi.M1[Xsi.M1$XsiParameter=="xsi.1",], Xsi.M2[Xsi.M2$XsiParameter=="xsi.1",], Xsi.M3[Xsi.M3$XsiParameter=="xsi.1",])

bp.trt <- ggplot(AOpts.df, aes(x=phys.grps, y=TrtRate, fill=Method)) + geom_boxplot() + labs(y="Optimal treatment rate", x="J") + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(text = element_text(size = 14))
bp.out <- ggplot(YOpts.df, aes(x=phys.grps, y=Outcome, fill=Method)) + geom_boxplot() + labs(y="Optimal outcome", x="J") + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(text = element_text(size = 14))
bp.psi0 <- ggplot(Psi.0.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\psi_{0}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi1 <- ggplot(Psi.1.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\psi_{1}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi2 <- ggplot(Psi.2.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\psi_{2}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=-1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi3 <- ggplot(Psi.3.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\psi_{3}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=-1.5, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.xsi0 <- ggplot(Xsi.0.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\xi_{0}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.xsi1 <- ggplot(Xsi.1.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\xi_{1}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))

ggsave(here::here("results", "sims_one_trt_and_trt_fr_mis_trtRate_bp.pdf"), plot = bp.trt)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_mis_outcome_bp.pdf"), plot = bp.out)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_mis_psi0_bp.pdf"), plot = bp.psi0)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_mis_psi1_bp.pdf"), plot = bp.psi1)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_mis_psi2_bp.pdf"), plot = bp.psi2)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_mis_psi3_bp.pdf"), plot = bp.psi3)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_mis_xsi0_bp.pdf"), plot = bp.xsi0)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_mis_xsi1_bp.pdf"), plot = bp.xsi1)

ggsave(here::here("results", "sims_one_trt_and_trt_fr_mis_trtRate_bp.jpg"), plot = bp.trt)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_mis_outcome_bp.jpg"), plot = bp.out)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_mis_psi0_bp.jpg"), plot = bp.psi0)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_mis_psi1_bp.jpg"), plot = bp.psi1)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_mis_psi2_bp.jpg"), plot = bp.psi2)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_mis_psi3_bp.jpg"), plot = bp.psi3)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_mis_xsi0_bp.jpg"), plot = bp.xsi0)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_mis_xsi1_bp.jpg"), plot = bp.xsi1)


#One Stage: Both Trt and Trt Free Models Correctly Specified

one_stage_trt_and_trt_fr_correct <- readRDS(here::here("results", "one_stage_trt_and_trt_fr_correct.rds"))

names(one_stage_trt_and_trt_fr_correct) <- c("AOpts_one_trt_and_trt_fr_correct", "YOpts_one_trt_and_trt_fr_correct", "Regret_one_trt_and_trt_fr_correct", "Psi_M1_one_trt_and_trt_fr_correct", "Psi_M2_one_trt_and_trt_fr_correct", "Psi_M3_one_trt_and_trt_fr_correct", "Xsi_M1_one_trt_and_trt_fr_correct", "Xsi_M2_one_trt_and_trt_fr_correct", "Xsi_M3_one_trt_and_trt_fr_correct", "Min_alters_one_trt_and_trt_fr_correct", "Max_alters_one_trt_and_trt_fr_correct", "Mean_alters_one_trt_and_trt_fr_correct", "Median_alters_one_trt_and_trt_fr_correct")

AOpts.df <- one_stage_trt_and_trt_fr_correct$AOpts_one_trt_and_trt_fr_correct
YOpts.df <- one_stage_trt_and_trt_fr_correct$YOpts_one_trt_and_trt_fr_correct
Regret.df <- one_stage_trt_and_trt_fr_correct$Regret_one_trt_and_trt_fr_correct
Psi.M1 <- one_stage_trt_and_trt_fr_correct$Psi_M1_one_trt_and_trt_fr_correct
Psi.M2 <- one_stage_trt_and_trt_fr_correct$Psi_M2_one_trt_and_trt_fr_correct
Psi.M3 <- one_stage_trt_and_trt_fr_correct$Psi_M3_one_trt_and_trt_fr_correct
Xsi.M1 <- one_stage_trt_and_trt_fr_correct$Xsi_M1_one_trt_and_trt_fr_correct
Xsi.M2 <- one_stage_trt_and_trt_fr_correct$Xsi_M2_one_trt_and_trt_fr_correct
Xsi.M3 <- one_stage_trt_and_trt_fr_correct$Xsi_M3_one_trt_and_trt_fr_correct

#Provides numerical, tabular summary of optimal treatment rates and regrets across methods
plot_trt_rates <- na.omit(AOpts.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(TrtRate),
         mean = mean(TrtRate),
         sd = sd(TrtRate),
         empse = sqrt((sum((TrtRate-mean(TrtRate))^2))/(n-1)),
         .groups="keep")

plot_trt_rates

plot_outcomes <- na.omit(YOpts.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(Outcome),
         mean = mean(Outcome),
         sd = sd(Outcome),
         empse = sqrt((sum((Outcome-mean(Outcome))^2))/(n-1)),
         .groups="keep")

plot_outcomes

plot_regrets <- na.omit(Regret.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(Regret),
         mean = mean(Regret),
         sd = sd(Regret),
         empse = sqrt((sum((Regret-mean(Regret))^2))/(n-1)),
         .groups="keep")

plot_regrets

write.csv(plot_trt_rates, here::here("results", "sims_one_trt_and_trt_fr_correct_trtRate_tbl.csv"))

write.csv(plot_outcomes, here::here("results", "sims_one_trt_and_trt_fr_correct_outcome_tbl.csv"))

write.csv(plot_regrets, here::here("results", "sims_one_trt_and_trt_fr_correct_regret_tbl.csv"))

AOpts.df$phys.grps <- as.factor(AOpts.df$phys.grps)
YOpts.df$phys.grps <- as.factor(YOpts.df$phys.grps)

Psi.M1$phys.grps <- as.factor(Psi.M1$phys.grps)
Psi.M1$Method <- rep("ME2", nrow(Psi.M1))

Psi.M2$phys.grps <- as.factor(Psi.M2$phys.grps)
Psi.M2$Method <- rep("ME3", nrow(Psi.M2))

Psi.M3$phys.grps <- as.factor(Psi.M3$phys.grps)
Psi.M3$Method <- rep("ME4", nrow(Psi.M3))

Xsi.M1$phys.grps <- as.factor(Xsi.M1$phys.grps)
Xsi.M1$Method <- rep("ME2", nrow(Xsi.M1))

Xsi.M2$phys.grps <- as.factor(Xsi.M2$phys.grps)
Xsi.M2$Method <- rep("ME3", nrow(Xsi.M2))

Xsi.M3$phys.grps <- as.factor(Xsi.M3$phys.grps)
Xsi.M3$Method <- rep("ME4", nrow(Xsi.M3))


Psi.0.df <- rbind(Psi.M1[Psi.M1$PsiParameter=="psi.0",], Psi.M2[Psi.M2$PsiParameter=="psi.0",], Psi.M3[Psi.M3$PsiParameter=="psi.0",])
Psi.1.df <- rbind(Psi.M1[Psi.M1$PsiParameter=="psi.1",], Psi.M2[Psi.M2$PsiParameter=="psi.1",], Psi.M3[Psi.M3$PsiParameter=="psi.1",])
Psi.2.df <- rbind(Psi.M1[Psi.M1$PsiParameter=="psi.2",], Psi.M2[Psi.M2$PsiParameter=="psi.2",], Psi.M3[Psi.M3$PsiParameter=="psi.2",])
Psi.3.df <- rbind(Psi.M1[Psi.M1$PsiParameter=="psi.3",], Psi.M2[Psi.M2$PsiParameter=="psi.3",], Psi.M3[Psi.M3$PsiParameter=="psi.3",])

Xsi.0.df <- rbind(Xsi.M1[Xsi.M1$XsiParameter=="xsi.0",], Xsi.M2[Xsi.M2$XsiParameter=="xsi.0",], Xsi.M3[Xsi.M3$XsiParameter=="xsi.0",])
Xsi.1.df <- rbind(Xsi.M1[Xsi.M1$XsiParameter=="xsi.1",], Xsi.M2[Xsi.M2$XsiParameter=="xsi.1",], Xsi.M3[Xsi.M3$XsiParameter=="xsi.1",])

bp.trt <- ggplot(AOpts.df, aes(x=phys.grps, y=TrtRate, fill=Method)) + geom_boxplot() + labs(y="Optimal treatment rate", x="J") + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(text = element_text(size = 14))
bp.out <- ggplot(YOpts.df, aes(x=phys.grps, y=Outcome, fill=Method)) + geom_boxplot() + labs(y="Optimal outcome", x="J") + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(text = element_text(size = 14))
bp.psi0 <- ggplot(Psi.0.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\psi_{0}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi1 <- ggplot(Psi.1.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\psi_{1}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi2 <- ggplot(Psi.2.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\psi_{2}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=-1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi3 <- ggplot(Psi.3.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\psi_{3}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=-1.5, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.xsi0 <- ggplot(Xsi.0.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\xi_{0}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.xsi1 <- ggplot(Xsi.1.df, aes(x=phys.grps, y=Estimate, fill=Method)) + geom_boxplot() + labs(x="J", y=TeX(r"($\xi_{1}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))

ggsave(here::here("results", "sims_one_trt_and_trt_fr_correct_trtRate_bp.pdf"), plot = bp.trt)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_correct_outcome_bp.pdf"), plot = bp.out)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_correct_psi0_bp.pdf"), plot = bp.psi0)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_correct_psi1_bp.pdf"), plot = bp.psi1)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_correct_psi2_bp.pdf"), plot = bp.psi2)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_correct_psi3_bp.pdf"), plot = bp.psi3)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_correct_xsi0_bp.pdf"), plot = bp.xsi0)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_correct_xsi1_bp.pdf"), plot = bp.xsi1)

ggsave(here::here("results", "sims_one_trt_and_trt_fr_correct_trtRate_bp.jpg"), plot = bp.trt)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_correct_outcome_bp.jpg"), plot = bp.out)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_correct_psi0_bp.jpg"), plot = bp.psi0)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_correct_psi1_bp.jpg"), plot = bp.psi1)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_correct_psi2_bp.jpg"), plot = bp.psi2)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_correct_psi3_bp.jpg"), plot = bp.psi3)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_correct_xsi0_bp.jpg"), plot = bp.xsi0)
ggsave(here::here("results", "sims_one_trt_and_trt_fr_correct_xsi1_bp.jpg"), plot = bp.xsi1)


#Two Stage: Trt Mis-specified

sims_two_stage_trt_mis <- readRDS(here::here("results", "sims_two_stage_trt_mis.rds"))

names(sims_two_stage_trt_mis) <- c("TwoStage_results_both_trt_mis", "TwoStage_results_both__resc_trt_mis", "TwoStage_results_wthn_trt_mis", "TwoStage_results_standard_trt_mis", "TwoStage_min_num_trt_mis", "TwoStage_max_num_trt_mis", "TwoStage_mean_num_trt_mis", "TwoStage_median_num_trt_mis")

run <- 1000
Results.both <- sims_two_stage_trt_mis$TwoStage_results_both_trt_mis
Results.both.resc <- sims_two_stage_trt_mis$TwoStage_results_both__resc_trt_mis
Results.wthn <- sims_two_stage_trt_mis$TwoStage_results_wthn_trt_mis
Results.standard <- sims_two_stage_trt_mis$TwoStage_results_standard_trt_mis

dwols.trt.df <- data.frame(AOpts = 100*c(c(Results.wthn$DWOLS.WTHN.STAGE1),c(Results.wthn$DWOLS.WTHN.STAGE2),c(Results.standard$DWOLS.STANDARD.STAGE1),c(Results.standard$DWOLS.STANDARD.STAGE2),c(Results.both$DWOLS.BOTH.STAGE1),c(Results.both$DWOLS.BOTH.STAGE2),c(Results.both.resc$DWOLS.BOTH.RESC.STAGE1),c(Results.both.resc$DWOLS.BOTH.RESC.STAGE2)), Method = rep(c("ME1","ME2","ME3", "ME4"),each=run*2), Stage = rep(rep(c("Stage 1","Stage 2"),each=run),4) )

bp.twostage.trt <- ggplot(dwols.trt.df, aes(x=Stage, y=AOpts, fill=Method)) + geom_boxplot() + xlab("Treatment Stage") + ylab("Optimal Treatment Rate (%)") + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(text = element_text(size = 14))
#Plot outcome results

dwols.out.df <- data.frame(YOpts = c(c(Results.both$Yopt), c(Results.wthn$Yopt.DWOLS.WTHN), c(Results.standard$Yopt.DWOLS.STANDARD), c(Results.both$Yopt.DWOLS.BOTH),c(Results.both.resc$Yopt.DWOLS.BOTH.RESC)), Outcome = rep(c("True", "ME1","ME2", "ME3","ME4"),each=1000))

dwols.out.df$Outcome <- factor(dwols.out.df$Outcome, levels=c("ME1","ME2", "ME3", "ME4","True"))
bp.twostage.out <- ggplot(dwols.out.df, aes(y=YOpts, fill=Outcome)) + 
  geom_boxplot() + scale_fill_grey(start=0.8, end=0.4) + ylab("Optimal Outcome") + theme(text = element_text(size = 14))

stage1.regret.df <- data.frame(Regret = c(c(Results.wthn$REGRET.STAGE1.WTHN), c(Results.standard$REGRET.STAGE1.STANDARD), c(Results.both$REGRET.STAGE1.BOTH),c(Results.both.resc$REGRET.STAGE1.BOTH.RESC)), Method = rep(c("ME1","ME2", "ME3","ME4"),each=1000))
stage2.regret.df <- data.frame(Regret = c(c(Results.wthn$REGRET.STAGE2.WTHN), c(Results.standard$REGRET.STAGE2.STANDARD), c(Results.both$REGRET.STAGE2.BOTH),c(Results.both.resc$REGRET.STAGE2.BOTH.RESC)), Method = rep(c("ME1","ME2", "ME3","ME4"),each=1000))
total.regret.df <- data.frame(Regret = c(c(Results.wthn$REGRET.TOTAL.WTHN), c(Results.standard$REGRET.TOTAL.STANDARD), c(Results.both$REGRET.TOTAL.BOTH),c(Results.both.resc$REGRET.TOTAL.BOTH.RESC)), Method = rep(c("ME1","ME2", "ME3","ME4"),each=1000))

plot_regrets_stage1 <- na.omit(stage1.regret.df) %>%
  group_by(Method) %>%
  summarize(n = n(),
         median = median(Regret),
         mean = mean(Regret),
         sd = sd(Regret),
         empse = sqrt((sum((Regret-mean(Regret))^2))/(n-1)),
         .groups="keep")

plot_regrets_stage1

plot_regrets_stage2 <- na.omit(stage2.regret.df) %>%
  group_by(Method) %>%
  summarize(n = n(),
         median = median(Regret),
         mean = mean(Regret),
         sd = sd(Regret),
         empse = sqrt((sum((Regret-mean(Regret))^2))/(n-1)),
         .groups="keep")

plot_regrets_stage2

plot_regrets_total <- na.omit(total.regret.df) %>%
  group_by(Method) %>%
  summarize(n = n(),
         median = median(Regret),
         mean = mean(Regret),
         sd = sd(Regret),
         empse = sqrt((sum((Regret-mean(Regret))^2))/(n-1)),
         .groups="keep")

plot_regrets_total

write.csv(plot_regrets_stage1, here::here("results", "sims_two_trt_mis_regret_stage1_tbl.csv"))
write.csv(plot_regrets_stage2, here::here("results", "sims_two_trt_mis_regret_stage2_tbl.csv"))
write.csv(plot_regrets_total, here::here("results", "sims_two_trt_mis_regret_total_tbl.csv"))

dwols.psi.stage1.df <- data.frame(Psi = c(c(Results.both$PSI.0.STAGE1.BOTH),c(Results.both$PSI.1.STAGE1.BOTH),c(Results.both$PSI.2.STAGE1.BOTH),c(Results.both$PSI.3.STAGE1.BOTH)), Blip = rep(c("Psi0", "Psi1", "Psi2", "Psi3"),each=1000))
dwols.psi.stage2.df <- data.frame(Psi = c(c(Results.both$PSI.0.STAGE2.BOTH),c(Results.both$PSI.1.STAGE2.BOTH),c(Results.both$PSI.2.STAGE2.BOTH),c(Results.both$PSI.3.STAGE2.BOTH)), Blip = rep(c("Psi0", "Psi1", "Psi2", "Psi3"),each=1000))

dwols.psi.stage1.df$Blip <- factor(dwols.psi.stage1.df$Blip, levels=c("Psi0", "Psi1", "Psi2", "Psi3"))
dwols.psi.stage2.df$Blip <- factor(dwols.psi.stage2.df$Blip, levels=c("Psi0", "Psi1", "Psi2", "Psi3"))
dwols.psi.stage1.df$TruePsi <- rep(c(1,1, -1, -1.5), each=1000)
dwols.psi.stage2.df$TruePsi <- rep(c(1,1, -1, -1.5), each=1000)

dwols.xsi.stage1.df <- data.frame(Xsi = c(c(Results.both$XSI.0.STAGE1.BOTH),c(Results.both$XSI.1.STAGE1.BOTH)), Spillover = rep(c("Xsi0", "Xsi1"),each=1000))
dwols.xsi.stage2.df <- data.frame(Xsi = c(c(Results.both$XSI.0.STAGE2.BOTH),c(Results.both$XSI.1.STAGE2.BOTH)), Spillover = rep(c("Xsi0", "Xsi1"),each=1000))

dwols.xsi.stage1.df$Spillover <- factor(dwols.xsi.stage1.df$Spillover, levels=c("Xsi0", "Xsi1"))
dwols.xsi.stage2.df$Spillover <- factor(dwols.xsi.stage2.df$Spillover, levels=c("Xsi0", "Xsi1"))
dwols.xsi.stage1.df$TrueXsi <- rep(c(1,1), each=1000)
dwols.xsi.stage2.df$TrueXsi <- rep(c(1,1), each=1000)

summary_xsi_dwols_both_1 <- na.omit(dwols.xsi.stage1.df) %>%
  group_by(Spillover) %>%
  summarize(n = n(),
         median = median(Xsi),
         mean = mean(Xsi),
         empse = sqrt((sum((Xsi-mean(Xsi))^2))/(n-1)),
         mse = (sum((Xsi-TrueXsi)^2))/n,
         .groups="keep")

summary_xsi_dwols_both_1

summary_xsi_dwols_both_2 <- na.omit(dwols.xsi.stage2.df) %>%
  group_by(Spillover) %>%
  summarize(n = n(),
         median = median(Xsi),
         mean = mean(Xsi),
         empse = sqrt((sum((Xsi-mean(Xsi))^2))/(n-1)),
         mse = (sum((Xsi-TrueXsi)^2))/n,
         .groups="keep")

summary_xsi_dwols_both_2

summary_psi_dwols_both_1 <- na.omit(dwols.psi.stage1.df) %>%
  group_by(Blip) %>%
  summarize(n = n(),
         median = median(Psi),
         mean = mean(Psi),
         empse = sqrt((sum((Psi-mean(Psi))^2))/(n-1)),
         mse = (sum((Psi-TruePsi)^2))/n,
         .groups="keep")

summary_psi_dwols_both_1

summary_psi_dwols_both_2 <- na.omit(dwols.psi.stage2.df) %>%
  group_by(Blip) %>%
  summarize(n = n(),
         median = median(Psi),
         mean = mean(Psi),
         empse = sqrt((sum((Psi-mean(Psi))^2))/(n-1)),
         mse = (sum((Psi-TruePsi)^2))/n,
         .groups="keep")

summary_psi_dwols_both_2

bp.stageone.psi <- ggplot(dwols.psi.stage1.df, aes(x=Blip, y=Psi, fill=Blip)) + 
  geom_boxplot(show.legend=FALSE) + xlab(TeX(r"($\psi$ Parameter)")) + ylab("Blip Parameters: Stage 1") + scale_x_discrete(labels=c("Psi0"=TeX(r"($\psi_{0}$)"), "Psi1"=TeX(r"($\psi_{1}$)"), "Psi2"=TeX(r"($\psi_{2}$)"), "Psi3"=TeX(r"($\psi_{3}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(aes(yintercept=1, linetype="1"), color = "red") + geom_hline(aes(yintercept=-1, linetype="-1"), color = "blue") + geom_hline(aes(yintercept=-1.5, linetype="-1.5"), color = "green")+ scale_linetype_manual(name ="True Value", values = c('dashed', 'dashed', 'dashed')) + theme(text = element_text(size = 14))

bp.stagetwo.psi <- ggplot(dwols.psi.stage2.df, aes(x=Blip, y=Psi, fill=Blip)) + 
  geom_boxplot(show.legend=FALSE) + xlab(TeX(r"($\psi$ Parameter)")) + ylab("Blip Parameters: Stage 2") + scale_x_discrete(labels=c("Psi0"=TeX(r"($\psi_{0}$)"), "Psi1"=TeX(r"($\psi_{1}$)"), "Psi2"=TeX(r"($\psi_{2}$)"), "Psi3"=TeX(r"($\psi_{3}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(aes(yintercept=1, linetype="1"), color = "red") + geom_hline(aes(yintercept=-1, linetype="-1"), color = "blue") + geom_hline(aes(yintercept=-1.5, linetype="-1.5"), color = "green")+ scale_linetype_manual(name ="True Value", values = c('dashed', 'dashed', 'dashed')) + theme(text = element_text(size = 14))

bp.stageone.xsi <- ggplot(dwols.xsi.stage1.df, aes(x=Spillover, y=Xsi, fill=Spillover)) + 
  geom_boxplot() + xlab(TeX(r"($\xi$ Parameter)")) + ylab("Spillover Parameters: Stage 1") + scale_x_discrete(labels=c("Xsi0"=TeX(r"($\xi_{0}$)"), "Xsi1"=TeX(r"($\xi_{1}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))

bp.stagetwo.xsi <- ggplot(dwols.xsi.stage2.df, aes(x=Spillover, y=Xsi, fill=Spillover)) + 
  geom_boxplot() + xlab(TeX(r"($\xi$ Parameter)")) + ylab("Spillover Parameters: Stage 2") + scale_x_discrete(labels=c("Xsi0"=TeX(r"($\xi_{0}$)"), "Xsi1"=TeX(r"($\xi_{1}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))


standard.psi.stage1.df <- data.frame(Psi = c(c(Results.standard$PSI.0.STAGE1.STANDARD),c(Results.standard$PSI.1.STAGE1.STANDARD),c(Results.standard$PSI.2.STAGE1.STANDARD),c(Results.standard$PSI.3.STAGE1.STANDARD)), Blip = rep(c("Psi0", "Psi1", "Psi2", "Psi3"),each=1000))
standard.psi.stage2.df <- data.frame(Psi = c(c(Results.standard$PSI.0.STAGE2.STANDARD),c(Results.standard$PSI.1.STAGE2.STANDARD),c(Results.standard$PSI.2.STAGE2.STANDARD),c(Results.standard$PSI.3.STAGE2.STANDARD)), Blip = rep(c("Psi0", "Psi1", "Psi2", "Psi3"),each=1000))

standard.psi.stage1.df$Blip <- factor(standard.psi.stage1.df$Blip, levels=c("Psi0", "Psi1", "Psi2", "Psi3"))
standard.psi.stage2.df$Blip <- factor(standard.psi.stage2.df$Blip, levels=c("Psi0", "Psi1", "Psi2", "Psi3"))
standard.psi.stage1.df$TruePsi <- rep(c(1, 1, -1, -1.5), each=1000)
standard.psi.stage2.df$TruePsi <- rep(c(1, 1, -1, -1.5), each=1000)

standard.xsi.stage1.df <- data.frame(Xsi = c(c(Results.standard$XSI.0.STAGE1.STANDARD),c(Results.standard$XSI.1.STAGE1.STANDARD)), Spillover = rep(c("Xsi0", "Xsi1"),each=1000))
standard.xsi.stage2.df <- data.frame(Xsi = c(c(Results.standard$XSI.0.STAGE2.STANDARD),c(Results.standard$XSI.1.STAGE2.STANDARD)), Spillover = rep(c("Xsi0", "Xsi1"),each=1000))

standard.xsi.stage1.df$Spillover <- factor(standard.xsi.stage1.df$Spillover, levels=c("Xsi0", "Xsi1"))
standard.xsi.stage2.df$Spillover <- factor(standard.xsi.stage2.df$Spillover, levels=c("Xsi0", "Xsi1"))
standard.xsi.stage1.df$TrueXsi <- rep(c(1,1), each=1000)
standard.xsi.stage2.df$TrueXsi <- rep(c(1,1), each=1000)

summary_xsi_standard_1 <- na.omit(standard.xsi.stage1.df) %>%
  group_by(Spillover) %>%
  summarize(n = n(),
         median = median(Xsi),
         mean = mean(Xsi),
         empse = sqrt((sum((Xsi-mean(Xsi))^2))/(n-1)),
         mse = (sum((Xsi-TrueXsi)^2))/n,
         .groups="keep")

summary_xsi_standard_1

summary_xsi_standard_2 <- na.omit(standard.xsi.stage2.df) %>%
  group_by(Spillover) %>%
  summarize(n = n(),
         median = median(Xsi),
         mean = mean(Xsi),
         empse = sqrt((sum((Xsi-mean(Xsi))^2))/(n-1)),
         mse = (sum((Xsi-TrueXsi)^2))/n,
         .groups="keep")

summary_xsi_standard_2

summary_psi_standard_1 <- na.omit(standard.psi.stage1.df) %>%
  group_by(Blip) %>%
  summarize(n = n(),
         median = median(Psi),
         mean = mean(Psi),
         empse = sqrt((sum((Psi-mean(Psi))^2))/(n-1)),
         mse = (sum((Psi-TruePsi)^2))/n,
         .groups="keep")

summary_psi_standard_1

summary_psi_standard_2 <- na.omit(standard.psi.stage2.df) %>%
  group_by(Blip) %>%
  summarize(n = n(),
         median = median(Psi),
         mean = mean(Psi),
         empse = sqrt((sum((Psi-mean(Psi))^2))/(n-1)),
         mse = (sum((Psi-TruePsi)^2))/n,
         .groups="keep")

summary_psi_standard_2


bp.stageone.psi.standard <- ggplot(standard.psi.stage1.df, aes(x=Blip, y=Psi, fill=Blip)) + 
  geom_boxplot(show.legend=FALSE) + xlab(TeX(r"($\psi$ Parameter)")) + ylab("Blip Parameters: Stage 1") + scale_x_discrete(labels=c("Psi0"=TeX(r"($\psi_{0}$)"), "Psi1"=TeX(r"($\psi_{1}$)"), "Psi2"=TeX(r"($\psi_{2}$)"), "Psi3"=TeX(r"($\psi_{3}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(aes(yintercept=1, linetype="1"), color = "red") + geom_hline(aes(yintercept=-1, linetype="-1"), color = "blue") + geom_hline(aes(yintercept=-1.5, linetype="-1.5"), color = "green")+ scale_linetype_manual(name ="True Value", values = c('dashed', 'dashed', 'dashed')) + theme(text = element_text(size = 14))

bp.stagetwo.psi.standard <- ggplot(standard.psi.stage2.df, aes(x=Blip, y=Psi, fill=Blip)) + 
  geom_boxplot(show.legend=FALSE) + xlab(TeX(r"($\psi$ Parameter)")) + ylab("Blip Parameters: Stage 2") + scale_x_discrete(labels=c("Psi0"=TeX(r"($\psi_{0}$)"), "Psi1"=TeX(r"($\psi_{1}$)"), "Psi2"=TeX(r"($\psi_{2}$)"), "Psi3"=TeX(r"($\psi_{3}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(aes(yintercept=1, linetype="1"), color = "red") + geom_hline(aes(yintercept=-1, linetype="-1"), color = "blue") + geom_hline(aes(yintercept=-1.5, linetype="-1.5"), color = "green")+ scale_linetype_manual(name ="True Value", values = c('dashed', 'dashed', 'dashed')) + theme(text = element_text(size = 14))

bp.stageone.xsi.standard <- ggplot(standard.xsi.stage1.df, aes(x=Spillover, y=Xsi, fill=Spillover)) + 
  geom_boxplot() + xlab(TeX(r"($\xi$ Parameter)")) + ylab("Spillover Parameters: Stage 1") + scale_x_discrete(labels=c("Xsi0"=TeX(r"($\xi_{0}$)"), "Xsi1"=TeX(r"($\xi_{1}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))

bp.stagetwo.xsi.standard <- ggplot(standard.xsi.stage2.df, aes(x=Spillover, y=Xsi, fill=Spillover)) + 
  geom_boxplot() + xlab(TeX(r"($\xi$ Parameter)")) + ylab("Spillover Parameters: Stage 2") + scale_x_discrete(labels=c("Xsi0"=TeX(r"($\xi_{0}$)"), "Xsi1"=TeX(r"($\xi_{1}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))

ggsave(here::here("results", "sims_two_trt_mis_trtRate_bp.pdf"), plot = bp.twostage.trt)
ggsave(here::here("results", "sims_two_trt_mis_outcome_bp.pdf"), plot = bp.twostage.out)
ggsave(here::here("results", "sims_two_trt_mis_stageone_psi_bp.pdf"), plot = bp.stageone.psi)
ggsave(here::here("results", "sims_two_trt_mis_stagetwo_psi_bp.pdf"), plot = bp.stagetwo.psi)
ggsave(here::here("results", "sims_two_trt_mis_stageone_xsi_bp.pdf"), plot = bp.stageone.xsi)
ggsave(here::here("results", "sims_two_trt_mis_stagetwo_xsi_bp.pdf"), plot = bp.stagetwo.xsi)
ggsave(here::here("results", "sims_two_trt_mis_stageone_psi_standard_bp.pdf"), plot = bp.stageone.psi.standard)
ggsave(here::here("results", "sims_two_trt_mis_stagetwo_psi_standard_bp.pdf"), plot = bp.stagetwo.psi.standard)
ggsave(here::here("results", "sims_two_trt_mis_stageone_xsi_standard_bp.pdf"), plot = bp.stageone.xsi.standard)
ggsave(here::here("results", "sims_two_trt_mis_stagetwo_xsi_standard_bp.pdf"), plot = bp.stagetwo.xsi.standard)

ggsave(here::here("results", "sims_two_trt_mis_trtRate_bp.jpg"), plot = bp.twostage.trt)
ggsave(here::here("results", "sims_two_trt_mis_outcome_bp.jpg"), plot = bp.twostage.out)
ggsave(here::here("results", "sims_two_trt_mis_stageone_psi_bp.jpg"), plot = bp.stageone.psi)
ggsave(here::here("results", "sims_two_trt_mis_stagetwo_psi_bp.jpg"), plot = bp.stagetwo.psi)
ggsave(here::here("results", "sims_two_trt_mis_stageone_xsi_bp.jpg"), plot = bp.stageone.xsi)
ggsave(here::here("results", "sims_two_trt_mis_stagetwo_xsi_bp.jpg"), plot = bp.stagetwo.xsi)
ggsave(here::here("results", "sims_two_trt_mis_stageone_psi_standard_bp.jpg"), plot = bp.stageone.psi.standard)
ggsave(here::here("results", "sims_two_trt_mis_stagetwo_psi_standard_bp.jpg"), plot = bp.stagetwo.psi.standard)
ggsave(here::here("results", "sims_two_trt_mis_stageone_xsi_standard_bp.jpg"), plot = bp.stageone.xsi.standard)
ggsave(here::here("results", "sims_two_trt_mis_stagetwo_xsi_standard_bp.jpg"), plot = bp.stagetwo.xsi.standard)


#Two Stage: Trt Free Mis-specified
sims_two_stage_trt_fr_mis <- readRDS(here::here("results", "sims_two_stage_trt_fr_mis.rds"))

names(sims_two_stage_trt_fr_mis) <- c("TwoStage_results_both_trt_fr_mis", "TwoStage_results_both__resc_trt_fr_mis", "TwoStage_results_wthn_trt_fr_mis", "TwoStage_results_standard_trt_fr_mis", "TwoStage_min_num_trt_fr_mis", "TwoStage_max_num_trt_fr_mis", "TwoStage_mean_num_trt_fr_mis", "TwoStage_median_num_trt_fr_mis")

run <- 1000
Results.both <- sims_two_stage_trt_fr_mis$TwoStage_results_both_trt_fr_mis
Results.both.resc <- sims_two_stage_trt_fr_mis$TwoStage_results_both__resc_trt_fr_mis
Results.wthn <- sims_two_stage_trt_fr_mis$TwoStage_results_wthn_trt_fr_mis
Results.standard <- sims_two_stage_trt_fr_mis$TwoStage_results_standard_trt_fr_mis

dwols.trt.df <- data.frame(AOpts = 100*c(c(Results.wthn$DWOLS.WTHN.STAGE1),c(Results.wthn$DWOLS.WTHN.STAGE2),c(Results.standard$DWOLS.STANDARD.STAGE1),c(Results.standard$DWOLS.STANDARD.STAGE2),c(Results.both$DWOLS.BOTH.STAGE1),c(Results.both$DWOLS.BOTH.STAGE2),c(Results.both.resc$DWOLS.BOTH.RESC.STAGE1),c(Results.both.resc$DWOLS.BOTH.RESC.STAGE2)), Method = rep(c("ME1","ME2","ME3", "ME4"),each=run*2), Stage = rep(rep(c("Stage 1","Stage 2"),each=run),4) )

bp.twostage.trt <- ggplot(dwols.trt.df, aes(x=Stage, y=AOpts, fill=Method)) + geom_boxplot() + xlab("Treatment Stage") + ylab("Optimal Treatment Rate (%)") + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(text = element_text(size = 14))
#Plot outcome results

dwols.out.df <- data.frame(YOpts = c(c(Results.both$Yopt), c(Results.wthn$Yopt.DWOLS.WTHN), c(Results.standard$Yopt.DWOLS.STANDARD), c(Results.both$Yopt.DWOLS.BOTH),c(Results.both.resc$Yopt.DWOLS.BOTH.RESC)), Outcome = rep(c("True", "ME1","ME2", "ME3","ME4"),each=1000))

dwols.out.df$Outcome <- factor(dwols.out.df$Outcome, levels=c("ME1","ME2", "ME3", "ME4","True"))
bp.twostage.out <- ggplot(dwols.out.df, aes(y=YOpts, fill=Outcome)) + 
  geom_boxplot() + scale_fill_grey(start=0.8, end=0.4) + ylab("Optimal Outcome") + theme(text = element_text(size = 14))

stage1.regret.df <- data.frame(Regret = c(c(Results.wthn$REGRET.STAGE1.WTHN), c(Results.standard$REGRET.STAGE1.STANDARD), c(Results.both$REGRET.STAGE1.BOTH),c(Results.both.resc$REGRET.STAGE1.BOTH.RESC)), Method = rep(c("ME1","ME2", "ME3","ME4"),each=1000))
stage2.regret.df <- data.frame(Regret = c(c(Results.wthn$REGRET.STAGE2.WTHN), c(Results.standard$REGRET.STAGE2.STANDARD), c(Results.both$REGRET.STAGE2.BOTH),c(Results.both.resc$REGRET.STAGE2.BOTH.RESC)), Method = rep(c("ME1","ME2", "ME3","ME4"),each=1000))
total.regret.df <- data.frame(Regret = c(c(Results.wthn$REGRET.TOTAL.WTHN), c(Results.standard$REGRET.TOTAL.STANDARD), c(Results.both$REGRET.TOTAL.BOTH),c(Results.both.resc$REGRET.TOTAL.BOTH.RESC)), Method = rep(c("ME1","ME2", "ME3","ME4"),each=1000))

plot_regrets_stage1 <- na.omit(stage1.regret.df) %>%
  group_by(Method) %>%
  summarize(n = n(),
         median = median(Regret),
         mean = mean(Regret),
         sd = sd(Regret),
         empse = sqrt((sum((Regret-mean(Regret))^2))/(n-1)),
         .groups="keep")

plot_regrets_stage1

plot_regrets_stage2 <- na.omit(stage2.regret.df) %>%
  group_by(Method) %>%
  summarize(n = n(),
         median = median(Regret),
         mean = mean(Regret),
         sd = sd(Regret),
         empse = sqrt((sum((Regret-mean(Regret))^2))/(n-1)),
         .groups="keep")

plot_regrets_stage2

plot_regrets_total <- na.omit(total.regret.df) %>%
  group_by(Method) %>%
  summarize(n = n(),
         median = median(Regret),
         mean = mean(Regret),
         sd = sd(Regret),
         empse = sqrt((sum((Regret-mean(Regret))^2))/(n-1)),
         .groups="keep")

plot_regrets_total

write.csv(plot_regrets_stage1, here::here("results", "sims_two_trt_fr_mis_regret_stage1_tbl.csv"))
write.csv(plot_regrets_stage2, here::here("results", "sims_two_trt_fr_mis_regret_stage2_tbl.csv"))
write.csv(plot_regrets_total, here::here("results", "sims_two_trt_fr_mis_regret_total_tbl.csv"))

dwols.psi.stage1.df <- data.frame(Psi = c(c(Results.both$PSI.0.STAGE1.BOTH),c(Results.both$PSI.1.STAGE1.BOTH),c(Results.both$PSI.2.STAGE1.BOTH),c(Results.both$PSI.3.STAGE1.BOTH)), Blip = rep(c("Psi0", "Psi1", "Psi2", "Psi3"),each=1000))
dwols.psi.stage2.df <- data.frame(Psi = c(c(Results.both$PSI.0.STAGE2.BOTH),c(Results.both$PSI.1.STAGE2.BOTH),c(Results.both$PSI.2.STAGE2.BOTH),c(Results.both$PSI.3.STAGE2.BOTH)), Blip = rep(c("Psi0", "Psi1", "Psi2", "Psi3"),each=1000))

dwols.psi.stage1.df$Blip <- factor(dwols.psi.stage1.df$Blip, levels=c("Psi0", "Psi1", "Psi2", "Psi3"))
dwols.psi.stage2.df$Blip <- factor(dwols.psi.stage2.df$Blip, levels=c("Psi0", "Psi1", "Psi2", "Psi3"))
dwols.psi.stage1.df$TruePsi <- rep(c(1,1, -1, -1.5), each=1000)
dwols.psi.stage2.df$TruePsi <- rep(c(1,1, -1, -1.5), each=1000)

dwols.xsi.stage1.df <- data.frame(Xsi = c(c(Results.both$XSI.0.STAGE1.BOTH),c(Results.both$XSI.1.STAGE1.BOTH)), Spillover = rep(c("Xsi0", "Xsi1"),each=1000))
dwols.xsi.stage2.df <- data.frame(Xsi = c(c(Results.both$XSI.0.STAGE2.BOTH),c(Results.both$XSI.1.STAGE2.BOTH)), Spillover = rep(c("Xsi0", "Xsi1"),each=1000))

dwols.xsi.stage1.df$Spillover <- factor(dwols.xsi.stage1.df$Spillover, levels=c("Xsi0", "Xsi1"))
dwols.xsi.stage2.df$Spillover <- factor(dwols.xsi.stage2.df$Spillover, levels=c("Xsi0", "Xsi1"))
dwols.xsi.stage1.df$TrueXsi <- rep(c(1,1), each=1000)
dwols.xsi.stage2.df$TrueXsi <- rep(c(1,1), each=1000)

summary_xsi_dwols_both_1 <- na.omit(dwols.xsi.stage1.df) %>%
  group_by(Spillover) %>%
  summarize(n = n(),
         median = median(Xsi),
         mean = mean(Xsi),
         empse = sqrt((sum((Xsi-mean(Xsi))^2))/(n-1)),
         mse = (sum((Xsi-TrueXsi)^2))/n,
         .groups="keep")

summary_xsi_dwols_both_1

summary_xsi_dwols_both_2 <- na.omit(dwols.xsi.stage2.df) %>%
  group_by(Spillover) %>%
  summarize(n = n(),
         median = median(Xsi),
         mean = mean(Xsi),
         empse = sqrt((sum((Xsi-mean(Xsi))^2))/(n-1)),
         mse = (sum((Xsi-TrueXsi)^2))/n,
         .groups="keep")

summary_xsi_dwols_both_2

summary_psi_dwols_both_1 <- na.omit(dwols.psi.stage1.df) %>%
  group_by(Blip) %>%
  summarize(n = n(),
         median = median(Psi),
         mean = mean(Psi),
         empse = sqrt((sum((Psi-mean(Psi))^2))/(n-1)),
         mse = (sum((Psi-TruePsi)^2))/n,
         .groups="keep")

summary_psi_dwols_both_1

summary_psi_dwols_both_2 <- na.omit(dwols.psi.stage2.df) %>%
  group_by(Blip) %>%
  summarize(n = n(),
         median = median(Psi),
         mean = mean(Psi),
         empse = sqrt((sum((Psi-mean(Psi))^2))/(n-1)),
         mse = (sum((Psi-TruePsi)^2))/n,
         .groups="keep")

summary_psi_dwols_both_2


bp.stageone.psi <- ggplot(dwols.psi.stage1.df, aes(x=Blip, y=Psi, fill=Blip)) + 
  geom_boxplot(show.legend=FALSE) + xlab(TeX(r"($\psi$ Parameter)")) + ylab("Blip Parameters: Stage 1") + scale_x_discrete(labels=c("Psi0"=TeX(r"($\psi_{0}$)"), "Psi1"=TeX(r"($\psi_{1}$)"), "Psi2"=TeX(r"($\psi_{2}$)"), "Psi3"=TeX(r"($\psi_{3}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(aes(yintercept=1, linetype="1"), color = "red") + geom_hline(aes(yintercept=-1, linetype="-1"), color = "blue") + geom_hline(aes(yintercept=-1.5, linetype="-1.5"), color = "green")+ scale_linetype_manual(name ="True Value", values = c('dashed', 'dashed', 'dashed')) + theme(text = element_text(size = 14))

bp.stagetwo.psi <- ggplot(dwols.psi.stage2.df, aes(x=Blip, y=Psi, fill=Blip)) + 
  geom_boxplot(show.legend=FALSE) + xlab(TeX(r"($\psi$ Parameter)")) + ylab("Blip Parameters: Stage 2") + scale_x_discrete(labels=c("Psi0"=TeX(r"($\psi_{0}$)"), "Psi1"=TeX(r"($\psi_{1}$)"), "Psi2"=TeX(r"($\psi_{2}$)"), "Psi3"=TeX(r"($\psi_{3}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(aes(yintercept=1, linetype="1"), color = "red") + geom_hline(aes(yintercept=-1, linetype="-1"), color = "blue") + geom_hline(aes(yintercept=-1.5, linetype="-1.5"), color = "green")+ scale_linetype_manual(name ="True Value", values = c('dashed', 'dashed', 'dashed')) + theme(text = element_text(size = 14))

bp.stageone.xsi <- ggplot(dwols.xsi.stage1.df, aes(x=Spillover, y=Xsi, fill=Spillover)) + 
  geom_boxplot() + xlab(TeX(r"($\xi$ Parameter)")) + ylab("Spillover Parameters: Stage 1") + scale_x_discrete(labels=c("Xsi0"=TeX(r"($\xi_{0}$)"), "Xsi1"=TeX(r"($\xi_{1}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))

bp.stagetwo.xsi <- ggplot(dwols.xsi.stage2.df, aes(x=Spillover, y=Xsi, fill=Spillover)) + 
  geom_boxplot() + xlab(TeX(r"($\xi$ Parameter)")) + ylab("Spillover Parameters: Stage 2") + scale_x_discrete(labels=c("Xsi0"=TeX(r"($\xi_{0}$)"), "Xsi1"=TeX(r"($\xi_{1}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))


standard.psi.stage1.df <- data.frame(Psi = c(c(Results.standard$PSI.0.STAGE1.STANDARD),c(Results.standard$PSI.1.STAGE1.STANDARD),c(Results.standard$PSI.2.STAGE1.STANDARD),c(Results.standard$PSI.3.STAGE1.STANDARD)), Blip = rep(c("Psi0", "Psi1", "Psi2", "Psi3"),each=1000))
standard.psi.stage2.df <- data.frame(Psi = c(c(Results.standard$PSI.0.STAGE2.STANDARD),c(Results.standard$PSI.1.STAGE2.STANDARD),c(Results.standard$PSI.2.STAGE2.STANDARD),c(Results.standard$PSI.3.STAGE2.STANDARD)), Blip = rep(c("Psi0", "Psi1", "Psi2", "Psi3"),each=1000))

standard.psi.stage1.df$Blip <- factor(standard.psi.stage1.df$Blip, levels=c("Psi0", "Psi1", "Psi2", "Psi3"))
standard.psi.stage2.df$Blip <- factor(standard.psi.stage2.df$Blip, levels=c("Psi0", "Psi1", "Psi2", "Psi3"))
standard.psi.stage1.df$TruePsi <- rep(c(1, 1, -1, -1.5), each=1000)
standard.psi.stage2.df$TruePsi <- rep(c(1, 1, -1, -1.5), each=1000)

standard.xsi.stage1.df <- data.frame(Xsi = c(c(Results.standard$XSI.0.STAGE1.STANDARD),c(Results.standard$XSI.1.STAGE1.STANDARD)), Spillover = rep(c("Xsi0", "Xsi1"),each=1000))
standard.xsi.stage2.df <- data.frame(Xsi = c(c(Results.standard$XSI.0.STAGE2.STANDARD),c(Results.standard$XSI.1.STAGE2.STANDARD)), Spillover = rep(c("Xsi0", "Xsi1"),each=1000))

standard.xsi.stage1.df$Spillover <- factor(standard.xsi.stage1.df$Spillover, levels=c("Xsi0", "Xsi1"))
standard.xsi.stage2.df$Spillover <- factor(standard.xsi.stage2.df$Spillover, levels=c("Xsi0", "Xsi1"))
standard.xsi.stage1.df$TrueXsi <- rep(c(1,1), each=1000)
standard.xsi.stage2.df$TrueXsi <- rep(c(1,1), each=1000)

summary_xsi_standard_1 <- na.omit(standard.xsi.stage1.df) %>%
  group_by(Spillover) %>%
  summarize(n = n(),
         median = median(Xsi),
         mean = mean(Xsi),
         empse = sqrt((sum((Xsi-mean(Xsi))^2))/(n-1)),
         mse = (sum((Xsi-TrueXsi)^2))/n,
         .groups="keep")

summary_xsi_standard_1

summary_xsi_standard_2 <- na.omit(standard.xsi.stage2.df) %>%
  group_by(Spillover) %>%
  summarize(n = n(),
         median = median(Xsi),
         mean = mean(Xsi),
         empse = sqrt((sum((Xsi-mean(Xsi))^2))/(n-1)),
         mse = (sum((Xsi-TrueXsi)^2))/n,
         .groups="keep")

summary_xsi_standard_2

summary_psi_standard_1 <- na.omit(standard.psi.stage1.df) %>%
  group_by(Blip) %>%
  summarize(n = n(),
         median = median(Psi),
         mean = mean(Psi),
         empse = sqrt((sum((Psi-mean(Psi))^2))/(n-1)),
         mse = (sum((Psi-TruePsi)^2))/n,
         .groups="keep")

summary_psi_standard_1

summary_psi_standard_2 <- na.omit(standard.psi.stage2.df) %>%
  group_by(Blip) %>%
  summarize(n = n(),
         median = median(Psi),
         mean = mean(Psi),
         empse = sqrt((sum((Psi-mean(Psi))^2))/(n-1)),
         mse = (sum((Psi-TruePsi)^2))/n,
         .groups="keep")

summary_psi_standard_2


bp.stageone.psi.standard <- ggplot(standard.psi.stage1.df, aes(x=Blip, y=Psi, fill=Blip)) + 
 geom_boxplot(show.legend=FALSE) + xlab(TeX(r"($\psi$ Parameter)")) + ylab("Blip Parameters: Stage 1") + scale_x_discrete(labels=c("Psi0"=TeX(r"($\psi_{0}$)"), "Psi1"=TeX(r"($\psi_{1}$)"), "Psi2"=TeX(r"($\psi_{2}$)"), "Psi3"=TeX(r"($\psi_{3}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(aes(yintercept=1, linetype="1"), color = "red") + geom_hline(aes(yintercept=-1, linetype="-1"), color = "blue") + geom_hline(aes(yintercept=-1.5, linetype="-1.5"), color = "green")+ scale_linetype_manual(name ="True Value", values = c('dashed', 'dashed', 'dashed')) + theme(text = element_text(size = 14))

bp.stagetwo.psi.standard <- ggplot(standard.psi.stage2.df, aes(x=Blip, y=Psi, fill=Blip)) + 
 geom_boxplot(show.legend=FALSE) + xlab(TeX(r"($\psi$ Parameter)")) + ylab("Blip Parameters: Stage 2") + scale_x_discrete(labels=c("Psi0"=TeX(r"($\psi_{0}$)"), "Psi1"=TeX(r"($\psi_{1}$)"), "Psi2"=TeX(r"($\psi_{2}$)"), "Psi3"=TeX(r"($\psi_{3}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(aes(yintercept=1, linetype="1"), color = "red") + geom_hline(aes(yintercept=-1, linetype="-1"), color = "blue") + geom_hline(aes(yintercept=-1.5, linetype="-1.5"), color = "green")+ scale_linetype_manual(name ="True Value", values = c('dashed', 'dashed', 'dashed')) + theme(text = element_text(size = 14))

bp.stageone.xsi.standard <- ggplot(standard.xsi.stage1.df, aes(x=Spillover, y=Xsi, fill=Spillover)) + 
  geom_boxplot() + xlab(TeX(r"($\xi$ Parameter)")) + ylab("Spillover Parameters: Stage 1") + scale_x_discrete(labels=c("Xsi0"=TeX(r"($\xi_{0}$)"), "Xsi1"=TeX(r"($\xi_{1}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))

bp.stagetwo.xsi.standard <- ggplot(standard.xsi.stage2.df, aes(x=Spillover, y=Xsi, fill=Spillover)) + 
  geom_boxplot() + xlab(TeX(r"($\xi$ Parameter)")) + ylab("Spillover Parameters: Stage 2") + scale_x_discrete(labels=c("Xsi0"=TeX(r"($\xi_{0}$)"), "Xsi1"=TeX(r"($\xi_{1}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))

ggsave(here::here("results", "sims_two_trt_fr_mis_trtRate_bp.pdf"), plot = bp.twostage.trt)
ggsave(here::here("results", "sims_two_trt_fr_mis_outcome_bp.pdf"), plot = bp.twostage.out)
ggsave(here::here("results", "sims_two_trt_fr_mis_stageone_psi_bp.pdf"), plot = bp.stageone.psi)
ggsave(here::here("results", "sims_two_trt_fr_mis_stagetwo_psi_bp.pdf"), plot = bp.stagetwo.psi)
ggsave(here::here("results", "sims_two_trt_fr_mis_stageone_xsi_bp.pdf"), plot = bp.stageone.xsi)
ggsave(here::here("results", "sims_two_trt_fr_mis_stagetwo_xsi_bp.pdf"), plot = bp.stagetwo.xsi)
ggsave(here::here("results", "sims_two_trt_fr_mis_stageone_psi_standard_bp.pdf"), plot = bp.stageone.psi.standard)
ggsave(here::here("results", "sims_two_trt_fr_mis_stagetwo_psi_standard_bp.pdf"), plot = bp.stagetwo.psi.standard)
ggsave(here::here("results", "sims_two_trt_fr_mis_stageone_xsi_standard_bp.pdf"), plot = bp.stageone.xsi.standard)
ggsave(here::here("results", "sims_two_trt_fr_mis_stagetwo_xsi_standard_bp.pdf"), plot = bp.stagetwo.xsi.standard)

ggsave(here::here("results", "sims_two_trt_fr_mis_trtRate_bp.jpg"), plot = bp.twostage.trt)
ggsave(here::here("results", "sims_two_trt_fr_mis_outcome_bp.jpg"), plot = bp.twostage.out)
ggsave(here::here("results", "sims_two_trt_fr_mis_stageone_psi_bp.jpg"), plot = bp.stageone.psi)
ggsave(here::here("results", "sims_two_trt_fr_mis_stagetwo_psi_bp.jpg"), plot = bp.stagetwo.psi)
ggsave(here::here("results", "sims_two_trt_fr_mis_stageone_xsi_bp.jpg"), plot = bp.stageone.xsi)
ggsave(here::here("results", "sims_two_trt_fr_mis_stagetwo_xsi_bp.jpg"), plot = bp.stagetwo.xsi)
ggsave(here::here("results", "sims_two_trt_fr_mis_stageone_psi_standard_bp.jpg"), plot = bp.stageone.psi.standard)
ggsave(here::here("results", "sims_two_trt_fr_mis_stagetwo_psi_standard_bp.jpg"), plot = bp.stagetwo.psi.standard)
ggsave(here::here("results", "sims_two_trt_fr_mis_stageone_xsi_standard_bp.jpg"), plot = bp.stageone.xsi.standard)
ggsave(here::here("results", "sims_two_trt_fr_mis_stagetwo_xsi_standard_bp.jpg"), plot = bp.stagetwo.xsi.standard)


#Two Stage: Trt Free Mis-specified Variation 2
sims_two_stage_trt_fr_mis_variation <- readRDS(here::here("results", "sims_two_stage_trt_fr_mis_variation.rds"))

names(sims_two_stage_trt_fr_mis_variation) <- c("TwoStage_results_both_trt_fr_mis_var", "TwoStage_results_both__resc_trt_fr_mis_var", "TwoStage_results_wthn_trt_fr_mis_var", "TwoStage_results_standard_trt_fr_mis_var", "TwoStage_min_num_trt_fr_mis_var", "TwoStage_max_num_trt_fr_mis_var", "TwoStage_mean_num_trt_fr_mis_var", "TwoStage_median_num_trt_fr_mis_var")

run <- 1000
Results.both <- sims_two_stage_trt_fr_mis_variation$TwoStage_results_both_trt_fr_mis_var
Results.both.resc <- sims_two_stage_trt_fr_mis_variation$TwoStage_results_both__resc_trt_fr_mis_var
Results.wthn <- sims_two_stage_trt_fr_mis_variation$TwoStage_results_wthn_trt_fr_mis_var
Results.standard <- sims_two_stage_trt_fr_mis_variation$TwoStage_results_standard_trt_fr_mis_var

dwols.trt.df <- data.frame(AOpts = 100*c(c(Results.wthn$DWOLS.WTHN.STAGE1),c(Results.wthn$DWOLS.WTHN.STAGE2),c(Results.standard$DWOLS.STANDARD.STAGE1),c(Results.standard$DWOLS.STANDARD.STAGE2),c(Results.both$DWOLS.BOTH.STAGE1),c(Results.both$DWOLS.BOTH.STAGE2),c(Results.both.resc$DWOLS.BOTH.RESC.STAGE1),c(Results.both.resc$DWOLS.BOTH.RESC.STAGE2)), Method = rep(c("ME1","ME2","ME3", "ME4"),each=run*2), Stage = rep(rep(c("Stage 1","Stage 2"),each=run),4) )

bp.twostage.trt <- ggplot(dwols.trt.df, aes(x=Stage, y=AOpts, fill=Method)) + geom_boxplot() + xlab("Treatment Stage") + ylab("Optimal Treatment Rate (%)") + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(text = element_text(size = 14))
#Plot outcome results

dwols.out.df <- data.frame(YOpts = c(c(Results.both$Yopt), c(Results.wthn$Yopt.DWOLS.WTHN), c(Results.standard$Yopt.DWOLS.STANDARD), c(Results.both$Yopt.DWOLS.BOTH),c(Results.both.resc$Yopt.DWOLS.BOTH.RESC)), Outcome = rep(c("True", "ME1","ME2", "ME3","ME4"),each=1000))

dwols.out.df$Outcome <- factor(dwols.out.df$Outcome, levels=c("ME1","ME2", "ME3", "ME4","True"))
bp.twostage.out <- ggplot(dwols.out.df, aes(y=YOpts, fill=Outcome)) + 
  geom_boxplot() + scale_fill_grey(start=0.8, end=0.4) + ylab("Optimal Outcome") + theme(text = element_text(size = 14))

stage1.regret.df <- data.frame(Regret = c(c(Results.wthn$REGRET.STAGE1.WTHN), c(Results.standard$REGRET.STAGE1.STANDARD), c(Results.both$REGRET.STAGE1.BOTH),c(Results.both.resc$REGRET.STAGE1.BOTH.RESC)), Method = rep(c("ME1","ME2", "ME3","ME4"),each=1000))
stage2.regret.df <- data.frame(Regret = c(c(Results.wthn$REGRET.STAGE2.WTHN), c(Results.standard$REGRET.STAGE2.STANDARD), c(Results.both$REGRET.STAGE2.BOTH),c(Results.both.resc$REGRET.STAGE2.BOTH.RESC)), Method = rep(c("ME1","ME2", "ME3","ME4"),each=1000))
total.regret.df <- data.frame(Regret = c(c(Results.wthn$REGRET.TOTAL.WTHN), c(Results.standard$REGRET.TOTAL.STANDARD), c(Results.both$REGRET.TOTAL.BOTH),c(Results.both.resc$REGRET.TOTAL.BOTH.RESC)), Method = rep(c("ME1","ME2", "ME3","ME4"),each=1000))

plot_regrets_stage1 <- na.omit(stage1.regret.df) %>%
  group_by(Method) %>%
  summarize(n = n(),
         median = median(Regret),
         mean = mean(Regret),
         sd = sd(Regret),
         empse = sqrt((sum((Regret-mean(Regret))^2))/(n-1)),
         .groups="keep")

plot_regrets_stage1

plot_regrets_stage2 <- na.omit(stage2.regret.df) %>%
  group_by(Method) %>%
  summarize(n = n(),
         median = median(Regret),
         mean = mean(Regret),
         sd = sd(Regret),
         empse = sqrt((sum((Regret-mean(Regret))^2))/(n-1)),
         .groups="keep")

plot_regrets_stage2

plot_regrets_total <- na.omit(total.regret.df) %>%
  group_by(Method) %>%
  summarize(n = n(),
         median = median(Regret),
         mean = mean(Regret),
         sd = sd(Regret),
         empse = sqrt((sum((Regret-mean(Regret))^2))/(n-1)),
         .groups="keep")

plot_regrets_total

write.csv(plot_regrets_stage1, here::here("results", "sims_two_trt_fr_mis_var_regret_stage1_tbl.csv"))
write.csv(plot_regrets_stage2, here::here("results", "sims_two_trt_fr_mis_var_regret_stage2_tbl.csv"))
write.csv(plot_regrets_total, here::here("results", "sims_two_trt_fr_mis_var_regret_total_tbl.csv"))


dwols.psi.stage1.df <- data.frame(Psi = c(c(Results.both$PSI.0.STAGE1.BOTH),c(Results.both$PSI.1.STAGE1.BOTH),c(Results.both$PSI.2.STAGE1.BOTH),c(Results.both$PSI.3.STAGE1.BOTH)), Blip = rep(c("Psi0", "Psi1", "Psi2", "Psi3"),each=1000))
dwols.psi.stage2.df <- data.frame(Psi = c(c(Results.both$PSI.0.STAGE2.BOTH),c(Results.both$PSI.1.STAGE2.BOTH),c(Results.both$PSI.2.STAGE2.BOTH),c(Results.both$PSI.3.STAGE2.BOTH)), Blip = rep(c("Psi0", "Psi1", "Psi2", "Psi3"),each=1000))

dwols.psi.stage1.df$Blip <- factor(dwols.psi.stage1.df$Blip, levels=c("Psi0", "Psi1", "Psi2", "Psi3"))
dwols.psi.stage2.df$Blip <- factor(dwols.psi.stage2.df$Blip, levels=c("Psi0", "Psi1", "Psi2", "Psi3"))
dwols.psi.stage1.df$TruePsi <- rep(c(1,1, -1, -1.5), each=1000)
dwols.psi.stage2.df$TruePsi <- rep(c(1,1, -1, -1.5), each=1000)

dwols.xsi.stage1.df <- data.frame(Xsi = c(c(Results.both$XSI.0.STAGE1.BOTH),c(Results.both$XSI.1.STAGE1.BOTH)), Spillover = rep(c("Xsi0", "Xsi1"),each=1000))
dwols.xsi.stage2.df <- data.frame(Xsi = c(c(Results.both$XSI.0.STAGE2.BOTH),c(Results.both$XSI.1.STAGE2.BOTH)), Spillover = rep(c("Xsi0", "Xsi1"),each=1000))

dwols.xsi.stage1.df$Spillover <- factor(dwols.xsi.stage1.df$Spillover, levels=c("Xsi0", "Xsi1"))
dwols.xsi.stage2.df$Spillover <- factor(dwols.xsi.stage2.df$Spillover, levels=c("Xsi0", "Xsi1"))
dwols.xsi.stage1.df$TrueXsi <- rep(c(1,1), each=1000)
dwols.xsi.stage2.df$TrueXsi <- rep(c(1,1), each=1000)

summary_xsi_dwols_both_1 <- na.omit(dwols.xsi.stage1.df) %>%
  group_by(Spillover) %>%
  summarize(n = n(),
         median = median(Xsi),
         mean = mean(Xsi),
         empse = sqrt((sum((Xsi-mean(Xsi))^2))/(n-1)),
         mse = (sum((Xsi-TrueXsi)^2))/n,
         .groups="keep")

summary_xsi_dwols_both_1

summary_xsi_dwols_both_2 <- na.omit(dwols.xsi.stage2.df) %>%
  group_by(Spillover) %>%
  summarize(n = n(),
         median = median(Xsi),
         mean = mean(Xsi),
         empse = sqrt((sum((Xsi-mean(Xsi))^2))/(n-1)),
         mse = (sum((Xsi-TrueXsi)^2))/n,
         .groups="keep")

summary_xsi_dwols_both_2

summary_psi_dwols_both_1 <- na.omit(dwols.psi.stage1.df) %>%
  group_by(Blip) %>%
  summarize(n = n(),
         median = median(Psi),
         mean = mean(Psi),
         empse = sqrt((sum((Psi-mean(Psi))^2))/(n-1)),
         mse = (sum((Psi-TruePsi)^2))/n,
         .groups="keep")

summary_psi_dwols_both_1

summary_psi_dwols_both_2 <- na.omit(dwols.psi.stage2.df) %>%
  group_by(Blip) %>%
  summarize(n = n(),
         median = median(Psi),
         mean = mean(Psi),
         empse = sqrt((sum((Psi-mean(Psi))^2))/(n-1)),
         mse = (sum((Psi-TruePsi)^2))/n,
         .groups="keep")

summary_psi_dwols_both_2


bp.stageone.psi <- ggplot(dwols.psi.stage1.df, aes(x=Blip, y=Psi, fill=Blip)) + 
  geom_boxplot(show.legend=FALSE) + xlab(TeX(r"($\psi$ Parameter)")) + ylab("Blip Parameters: Stage 1") + scale_x_discrete(labels=c("Psi0"=TeX(r"($\psi_{0}$)"), "Psi1"=TeX(r"($\psi_{1}$)"), "Psi2"=TeX(r"($\psi_{2}$)"), "Psi3"=TeX(r"($\psi_{3}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(aes(yintercept=1, linetype="1"), color = "red") + geom_hline(aes(yintercept=-1, linetype="-1"), color = "blue") + geom_hline(aes(yintercept=-1.5, linetype="-1.5"), color = "green")+ scale_linetype_manual(name ="True Value", values = c('dashed', 'dashed', 'dashed')) + theme(text = element_text(size = 14))

bp.stagetwo.psi <-ggplot(dwols.psi.stage2.df, aes(x=Blip, y=Psi, fill=Blip)) + 
  geom_boxplot(show.legend=FALSE) + xlab(TeX(r"($\psi$ Parameter)")) + ylab("Blip Parameters: Stage 2") + scale_x_discrete(labels=c("Psi0"=TeX(r"($\psi_{0}$)"), "Psi1"=TeX(r"($\psi_{1}$)"), "Psi2"=TeX(r"($\psi_{2}$)"), "Psi3"=TeX(r"($\psi_{3}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(aes(yintercept=1, linetype="1"), color = "red") + geom_hline(aes(yintercept=-1, linetype="-1"), color = "blue") + geom_hline(aes(yintercept=-1.5, linetype="-1.5"), color = "green")+ scale_linetype_manual(name ="True Value", values = c('dashed', 'dashed', 'dashed')) + theme(text = element_text(size = 14))

bp.stageone.xsi <-ggplot(dwols.xsi.stage1.df, aes(x=Spillover, y=Xsi, fill=Spillover)) + 
  geom_boxplot() + xlab(TeX(r"($\xi$ Parameter)")) + ylab("Spillover Parameters: Stage 1") + scale_x_discrete(labels=c("Xsi0"=TeX(r"($\xi_{0}$)"), "Xsi1"=TeX(r"($\xi_{1}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))

bp.stagetwo.xsi <-ggplot(dwols.xsi.stage2.df, aes(x=Spillover, y=Xsi, fill=Spillover)) + 
  geom_boxplot() + xlab(TeX(r"($\xi$ Parameter)")) + ylab("Spillover Parameters: Stage 2") + scale_x_discrete(labels=c("Xsi0"=TeX(r"($\xi_{0}$)"), "Xsi1"=TeX(r"($\xi_{1}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))


standard.psi.stage1.df <- data.frame(Psi = c(c(Results.standard$PSI.0.STAGE1.STANDARD),c(Results.standard$PSI.1.STAGE1.STANDARD),c(Results.standard$PSI.2.STAGE1.STANDARD),c(Results.standard$PSI.3.STAGE1.STANDARD)), Blip = rep(c("Psi0", "Psi1", "Psi2", "Psi3"),each=1000))
standard.psi.stage2.df <- data.frame(Psi = c(c(Results.standard$PSI.0.STAGE2.STANDARD),c(Results.standard$PSI.1.STAGE2.STANDARD),c(Results.standard$PSI.2.STAGE2.STANDARD),c(Results.standard$PSI.3.STAGE2.STANDARD)), Blip = rep(c("Psi0", "Psi1", "Psi2", "Psi3"),each=1000))

standard.psi.stage1.df$Blip <- factor(standard.psi.stage1.df$Blip, levels=c("Psi0", "Psi1", "Psi2", "Psi3"))
standard.psi.stage2.df$Blip <- factor(standard.psi.stage2.df$Blip, levels=c("Psi0", "Psi1", "Psi2", "Psi3"))
standard.psi.stage1.df$TruePsi <- rep(c(1, 1, -1, -1.5), each=1000)
standard.psi.stage2.df$TruePsi <- rep(c(1, 1, -1, -1.5), each=1000)

standard.xsi.stage1.df <- data.frame(Xsi = c(c(Results.standard$XSI.0.STAGE1.STANDARD),c(Results.standard$XSI.1.STAGE1.STANDARD)), Spillover = rep(c("Xsi0", "Xsi1"),each=1000))
standard.xsi.stage2.df <- data.frame(Xsi = c(c(Results.standard$XSI.0.STAGE2.STANDARD),c(Results.standard$XSI.1.STAGE2.STANDARD)), Spillover = rep(c("Xsi0", "Xsi1"),each=1000))

standard.xsi.stage1.df$Spillover <- factor(standard.xsi.stage1.df$Spillover, levels=c("Xsi0", "Xsi1"))
standard.xsi.stage2.df$Spillover <- factor(standard.xsi.stage2.df$Spillover, levels=c("Xsi0", "Xsi1"))
standard.xsi.stage1.df$TrueXsi <- rep(c(1,1), each=1000)
standard.xsi.stage2.df$TrueXsi <- rep(c(1,1), each=1000)

summary_xsi_standard_1 <- na.omit(standard.xsi.stage1.df) %>%
  group_by(Spillover) %>%
  summarize(n = n(),
         median = median(Xsi),
         mean = mean(Xsi),
         empse = sqrt((sum((Xsi-mean(Xsi))^2))/(n-1)),
         mse = (sum((Xsi-TrueXsi)^2))/n,
         .groups="keep")

summary_xsi_standard_1

summary_xsi_standard_2 <- na.omit(standard.xsi.stage2.df) %>%
  group_by(Spillover) %>%
  summarize(n = n(),
         median = median(Xsi),
         mean = mean(Xsi),
         empse = sqrt((sum((Xsi-mean(Xsi))^2))/(n-1)),
         mse = (sum((Xsi-TrueXsi)^2))/n,
         .groups="keep")

summary_xsi_standard_2

summary_psi_standard_1 <- na.omit(standard.psi.stage1.df) %>%
  group_by(Blip) %>%
  summarize(n = n(),
         median = median(Psi),
         mean = mean(Psi),
         empse = sqrt((sum((Psi-mean(Psi))^2))/(n-1)),
         mse = (sum((Psi-TruePsi)^2))/n,
         .groups="keep")

summary_psi_standard_1

summary_psi_standard_2 <- na.omit(standard.psi.stage2.df) %>%
  group_by(Blip) %>%
  summarize(n = n(),
         median = median(Psi),
         mean = mean(Psi),
         empse = sqrt((sum((Psi-mean(Psi))^2))/(n-1)),
         mse = (sum((Psi-TruePsi)^2))/n,
         .groups="keep")

summary_psi_standard_2


bp.stageone.psi.standard <-ggplot(standard.psi.stage1.df, aes(x=Blip, y=Psi, fill=Blip)) + 
  geom_boxplot(show.legend=FALSE) + xlab(TeX(r"($\psi$ Parameter)")) + ylab("Blip Parameters: Stage 1") + scale_x_discrete(labels=c("Psi0"=TeX(r"($\psi_{0}$)"), "Psi1"=TeX(r"($\psi_{1}$)"), "Psi2"=TeX(r"($\psi_{2}$)"), "Psi3"=TeX(r"($\psi_{3}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(aes(yintercept=1, linetype="1"), color = "red") + geom_hline(aes(yintercept=-1, linetype="-1"), color = "blue") + geom_hline(aes(yintercept=-1.5, linetype="-1.5"), color = "green")+ scale_linetype_manual(name ="True Value", values = c('dashed', 'dashed', 'dashed')) + theme(text = element_text(size = 14))

bp.stagetwo.psi.standard <-ggplot(standard.psi.stage2.df, aes(x=Blip, y=Psi, fill=Blip)) + 
  geom_boxplot(show.legend=FALSE) + xlab(TeX(r"($\psi$ Parameter)")) + ylab("Blip Parameters: Stage 2") + scale_x_discrete(labels=c("Psi0"=TeX(r"($\psi_{0}$)"), "Psi1"=TeX(r"($\psi_{1}$)"), "Psi2"=TeX(r"($\psi_{2}$)"), "Psi3"=TeX(r"($\psi_{3}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + geom_hline(aes(yintercept=1, linetype="1"), color = "red") + geom_hline(aes(yintercept=-1, linetype="-1"), color = "blue") + geom_hline(aes(yintercept=-1.5, linetype="-1.5"), color = "green")+ scale_linetype_manual(name ="True Value", values = c('dashed', 'dashed', 'dashed')) + theme(text = element_text(size = 14))

bp.stageone.xsi.standard <-ggplot(standard.xsi.stage1.df, aes(x=Spillover, y=Xsi, fill=Spillover)) + 
  geom_boxplot() + xlab(TeX(r"($\xi$ Parameter)")) + ylab("Spillover Parameters: Stage 1") + scale_x_discrete(labels=c("Xsi0"=TeX(r"($\xi_{0}$)"), "Xsi1"=TeX(r"($\xi_{1}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))

bp.stagetwo.xsi.standard <-ggplot(standard.xsi.stage2.df, aes(x=Spillover, y=Xsi, fill=Spillover)) + 
  geom_boxplot() + xlab(TeX(r"($\xi$ Parameter)")) + ylab("Spillover Parameters: Stage 2") + scale_x_discrete(labels=c("Xsi0"=TeX(r"($\xi_{0}$)"), "Xsi1"=TeX(r"($\xi_{1}$)"))) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))

ggsave(here::here("results", "sims_two_trt_fr_mis_var_trtRate_bp.pdf"), plot = bp.twostage.trt)
ggsave(here::here("results", "sims_two_trt_fr_mis_var_outcome_bp.pdf"), plot = bp.twostage.out)
ggsave(here::here("results", "sims_two_trt_fr_mis_var_stageone_psi_bp.pdf"), plot = bp.stageone.psi)
ggsave(here::here("results", "sims_two_trt_fr_mis_var_stagetwo_psi_bp.pdf"), plot = bp.stagetwo.psi)
ggsave(here::here("results", "sims_two_trt_fr_mis_var_stageone_xsi_bp.pdf"), plot = bp.stageone.xsi)
ggsave(here::here("results", "sims_two_trt_fr_mis_var_stagetwo_xsi_bp.pdf"), plot = bp.stagetwo.xsi)
ggsave(here::here("results", "sims_two_trt_fr_mis_var_stageone_psi_standard_bp.pdf"), plot = bp.stageone.psi.standard)
ggsave(here::here("results", "sims_two_trt_fr_mis_var_stagetwo_psi_standard_bp.pdf"), plot = bp.stagetwo.psi.standard)
ggsave(here::here("results", "sims_two_trt_fr_mis_var_stageone_xsi_standard_bp.pdf"), plot = bp.stageone.xsi.standard)
ggsave(here::here("results", "sims_two_trt_fr_mis_var_stagetwo_xsi_standard_bp.pdf"), plot = bp.stagetwo.xsi.standard)

ggsave(here::here("results", "sims_two_trt_fr_mis_var_trtRate_bp.jpg"), plot = bp.twostage.trt)
ggsave(here::here("results", "sims_two_trt_fr_mis_var_outcome_bp.jpg"), plot = bp.twostage.out)
ggsave(here::here("results", "sims_two_trt_fr_mis_var_stageone_psi_bp.jpg"), plot = bp.stageone.psi)
ggsave(here::here("results", "sims_two_trt_fr_mis_var_stagetwo_psi_bp.jpg"), plot = bp.stagetwo.psi)
ggsave(here::here("results", "sims_two_trt_fr_mis_var_stageone_xsi_bp.jpg"), plot = bp.stageone.xsi)
ggsave(here::here("results", "sims_two_trt_fr_mis_var_stagetwo_xsi_bp.jpg"), plot = bp.stagetwo.xsi)
ggsave(here::here("results", "sims_two_trt_fr_mis_var_stageone_psi_standard_bp.jpg"), plot = bp.stageone.psi.standard)
ggsave(here::here("results", "sims_two_trt_fr_mis_var_stagetwo_psi_standard_bp.jpg"), plot = bp.stagetwo.psi.standard)
ggsave(here::here("results", "sims_two_trt_fr_mis_var_stageone_xsi_standard_bp.jpg"), plot = bp.stageone.xsi.standard)
ggsave(here::here("results", "sims_two_trt_fr_mis_var_stagetwo_xsi_standard_bp.jpg"), plot = bp.stagetwo.xsi.standard)


#***************************************************************************************************


#One Stage: ego-alter treatment-free model components mis-specified, same with treatment model

study1b_extension_trt_mis <- readRDS(here::here("results", "study1b_extension_trt_mis.rds"))

names(study1b_extension_trt_mis) <- c("AOpts_study1b_ext_trt_mis", "YOpts_study1b_ext_trt_mis", "Regret_study1b_ext_trt_mis", "Psi_Dw_study1b_ext_trt_mis", "Psi_M1_study1b_ext_trt_mis", "Psi_M2_study1b_ext_trt_mis", "Psi_M3_study1b_ext_trt_mis", "Psi_M4_study1b_ext_trt_mis", "Psi_M5_study1b_ext_trt_mis", "Xsi_Dw_study1b_ext_trt_mis", "Xsi_M1_study1b_ext_trt_mis", "Xsi_M2_study1b_ext_trt_mis", "Xsi_M3_study1b_ext_trt_mis", "Xsi_M4_study1b_ext_trt_mis", "Xsi_M5_study1b_ext_trt_mis", "Mean_W_alters_study1b_ext_trt_mis", "Mean_B_alters_study1b_ext_trt_mis")

AOpts.df <- study1b_extension_trt_mis$AOpts_study1b_ext_trt_mis
YOpts.df <- study1b_extension_trt_mis$YOpts_study1b_ext_trt_mis
Regret.df <- study1b_extension_trt_mis$Regret_study1b_ext_trt_mis

Psi.Dw <- study1b_extension_trt_mis$Psi_Dw_study1b_ext_trt_mis
Psi.M1 <- study1b_extension_trt_mis$Psi_M1_study1b_ext_trt_mis
Psi.M2 <- study1b_extension_trt_mis$Psi_M2_study1b_ext_trt_mis
Psi.M3 <- study1b_extension_trt_mis$Psi_M3_study1b_ext_trt_mis
Psi.M4 <- study1b_extension_trt_mis$Psi_M4_study1b_ext_trt_mis
Psi.M5 <- study1b_extension_trt_mis$Psi_M5_study1b_ext_trt_mis

Xsi.Dw <- study1b_extension_trt_mis$Xsi_Dw_study1b_ext_trt_mis
Xsi.M1 <- study1b_extension_trt_mis$Xsi_M1_study1b_ext_trt_mis
Xsi.M2 <- study1b_extension_trt_mis$Xsi_M2_study1b_ext_trt_mis
Xsi.M3 <- study1b_extension_trt_mis$Xsi_M3_study1b_ext_trt_mis
Xsi.M4 <- study1b_extension_trt_mis$Xsi_M4_study1b_ext_trt_mis
Xsi.M5 <- study1b_extension_trt_mis$Xsi_M5_study1b_ext_trt_mis

#Provides numerical, tabular summary of optimal treatment rates and regrets across methods
plot_trt_rates <- na.omit(AOpts.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(TrtRate),
         mean = mean(TrtRate),
         sd = sd(TrtRate),
         empse = sqrt((sum((TrtRate-mean(TrtRate))^2))/(n-1)),
         .groups="keep")

plot_trt_rates

plot_outcomes <- na.omit(YOpts.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(Outcome),
         mean = mean(Outcome),
         sd = sd(Outcome),
         empse = sqrt((sum((Outcome-mean(Outcome))^2))/(n-1)),
         .groups="keep")

plot_outcomes

plot_regretss <- na.omit(Regret.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(Regret),
         mean = mean(Regret),
         sd = sd(Regret),
         empse = sqrt((sum((Regret-mean(Regret))^2))/(n-1)),
         .groups="keep")

plot_regrets

write.csv(plot_trt_rates, here::here("results", "sims_1b_extension_trt_mis_trtRate_tbl.csv"))

write.csv(plot_outcomes, here::here("results", "sims_1b_extension_trt_mis_outcome_tbl.csv"))

write.csv(plot_regrets, here::here("results", "sims_1b_extension_trt_mis_regret_tbl.csv"))

AOpts.df$phys.grps <- as.factor(AOpts.df$phys.grps)
YOpts.df$phys.grps <- as.factor(YOpts.df$phys.grps)

Psi.Dw$phys.grps <- as.factor(Psi.Dw$phys.grps)
Psi.Dw$Method <- rep("ME1", nrow(Psi.Dw))

Psi.M1$phys.grps <- as.factor(Psi.M1$phys.grps)
Psi.M1$Method <- rep("ME2", nrow(Psi.M1))

Psi.M2$phys.grps <- as.factor(Psi.M2$phys.grps)
Psi.M2$Method <- rep("ME3", nrow(Psi.M2))

Psi.M3$phys.grps <- as.factor(Psi.M3$phys.grps)
Psi.M3$Method <- rep("ME4", nrow(Psi.M3))

Psi.M4$phys.grps <- as.factor(Psi.M4$phys.grps)
Psi.M4$Method <- rep("ME5", nrow(Psi.M4))

Psi.M5$phys.grps <- as.factor(Psi.M5$phys.grps)
Psi.M5$Method <- rep("ME6", nrow(Psi.M5))

Xsi.Dw$phys.grps <- as.factor(Xsi.Dw$phys.grps)
Xsi.Dw$Method <- rep("ME1", nrow(Xsi.Dw))

Xsi.M1$phys.grps <- as.factor(Xsi.M1$phys.grps)
Xsi.M1$Method <- rep("ME2", nrow(Xsi.M1))

Xsi.M2$phys.grps <- as.factor(Xsi.M2$phys.grps)
Xsi.M2$Method <- rep("ME3", nrow(Xsi.M2))

Xsi.M3$phys.grps <- as.factor(Xsi.M3$phys.grps)
Xsi.M3$Method <- rep("ME4", nrow(Xsi.M3))

Xsi.M4$phys.grps <- as.factor(Xsi.M4$phys.grps)
Xsi.M4$Method <- rep("ME5", nrow(Xsi.M4))

Xsi.M5$phys.grps <- as.factor(Xsi.M5$phys.grps)
Xsi.M5$Method <- rep("ME6", nrow(Xsi.M5))


Psi.0.df <- rbind(Psi.Dw[Psi.Dw$PsiParameter=="psi.0",], Psi.M1[Psi.M1$PsiParameter=="psi.0",], Psi.M2[Psi.M2$PsiParameter=="psi.0",], Psi.M3[Psi.M3$PsiParameter=="psi.0",], Psi.M4[Psi.M4$PsiParameter=="psi.0",], Psi.M5[Psi.M5$PsiParameter=="psi.0",])
Psi.1.df <- rbind(Psi.Dw[Psi.Dw$PsiParameter=="psi.1",], Psi.M1[Psi.M1$PsiParameter=="psi.1",], Psi.M2[Psi.M2$PsiParameter=="psi.1",], Psi.M3[Psi.M3$PsiParameter=="psi.1",], Psi.M4[Psi.M4$PsiParameter=="psi.1",], Psi.M5[Psi.M5$PsiParameter=="psi.1",])
Psi.2.df <- rbind(Psi.Dw[Psi.Dw$PsiParameter=="psi.2",], Psi.M1[Psi.M1$PsiParameter=="psi.2",], Psi.M2[Psi.M2$PsiParameter=="psi.2",], Psi.M3[Psi.M3$PsiParameter=="psi.2",], Psi.M4[Psi.M4$PsiParameter=="psi.2",], Psi.M5[Psi.M5$PsiParameter=="psi.2",])
Psi.3.df <- rbind(Psi.Dw[Psi.Dw$PsiParameter=="psi.3",], Psi.M1[Psi.M1$PsiParameter=="psi.3",], Psi.M2[Psi.M2$PsiParameter=="psi.3",], Psi.M3[Psi.M3$PsiParameter=="psi.3",], Psi.M4[Psi.M4$PsiParameter=="psi.3",], Psi.M5[Psi.M5$PsiParameter=="psi.3",])

Xsi.0.df <- rbind(Xsi.Dw[Xsi.Dw$XsiParameter=="xsi.0",], Xsi.M1[Xsi.M1$XsiParameter=="xsi.0",], Xsi.M2[Xsi.M2$XsiParameter=="xsi.0",], Xsi.M3[Xsi.M3$XsiParameter=="xsi.0",], Xsi.M4[Xsi.M4$XsiParameter=="xsi.0",], Xsi.M5[Xsi.M5$XsiParameter=="xsi.0",])
Xsi.1.df <- rbind(Xsi.Dw[Xsi.Dw$XsiParameter=="xsi.1",], Xsi.M1[Xsi.M1$XsiParameter=="xsi.1",], Xsi.M2[Xsi.M2$XsiParameter=="xsi.1",], Xsi.M3[Xsi.M3$XsiParameter=="xsi.1",], Xsi.M4[Xsi.M4$XsiParameter=="xsi.1",], Xsi.M5[Xsi.M5$XsiParameter=="xsi.1",])

bp.trt <- ggplot(AOpts.df, aes(x=Method, y=TrtRate, fill=Method)) + geom_boxplot() + labs(y="Optimal treatment rate") + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + theme(text = element_text(size = 14))
bp.out <- ggplot(YOpts.df, aes(x=Method, y=Outcome, fill=Method)) + geom_boxplot() + labs(y="Optimal outcome") + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + theme(text = element_text(size = 14))
bp.psi0 <- ggplot(Psi.0.df, aes(x=Method, y=Estimate, fill=Method)) + geom_boxplot() + labs(y=TeX(r"($\psi_{0}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi1 <- ggplot(Psi.1.df, aes(x=Method, y=Estimate, fill=Method)) + geom_boxplot() + labs(y=TeX(r"($\psi_{1}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi2 <- ggplot(Psi.2.df, aes(x=Method, y=Estimate, fill=Method)) + geom_boxplot() + labs(y=TeX(r"($\psi_{2}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=-1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi3 <- ggplot(Psi.3.df, aes(x=Method, y=Estimate, fill=Method)) + geom_boxplot() + labs(y=TeX(r"($\psi_{3}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=-1.5, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.xsi0 <- ggplot(Xsi.0.df, aes(x=Method, y=Estimate, fill=Method)) + geom_boxplot() + labs(y=TeX(r"($\xi_{0}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.xsi1 <- ggplot(Xsi.1.df, aes(x=Method, y=Estimate, fill=Method)) + geom_boxplot() + labs(y=TeX(r"($\xi_{1}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))

ggsave(here::here("results", "sims_1b_extension_trt_mis_trtRate_bp.pdf"), plot = bp.trt)
ggsave(here::here("results", "sims_1b_extension_trt_mis_outcome_bp.pdf"), plot = bp.out)
ggsave(here::here("results", "sims_1b_extension_trt_mis_psi0_bp.pdf"), plot = bp.psi0)
ggsave(here::here("results", "sims_1b_extension_trt_mis_psi1_bp.pdf"), plot = bp.psi1)
ggsave(here::here("results", "sims_1b_extension_trt_mis_psi2_bp.pdf"), plot = bp.psi2)
ggsave(here::here("results", "sims_1b_extension_trt_mis_psi3_bp.pdf"), plot = bp.psi3)
ggsave(here::here("results", "sims_1b_extension_trt_mis_xsi0_bp.pdf"), plot = bp.xsi0)
ggsave(here::here("results", "sims_1b_extension_trt_mis_xsi1_bp.pdf"), plot = bp.xsi1)

ggsave(here::here("results", "sims_1b_extension_trt_mis_trtRate_bp.jpg"), plot = bp.trt)
ggsave(here::here("results", "sims_1b_extension_trt_mis_outcome_bp.jpg"), plot = bp.out)
ggsave(here::here("results", "sims_1b_extension_trt_mis_psi0_bp.jpg"), plot = bp.psi0)
ggsave(here::here("results", "sims_1b_extension_trt_mis_psi1_bp.jpg"), plot = bp.psi1)
ggsave(here::here("results", "sims_1b_extension_trt_mis_psi2_bp.jpg"), plot = bp.psi2)
ggsave(here::here("results", "sims_1b_extension_trt_mis_psi3_bp.jpg"), plot = bp.psi3)
ggsave(here::here("results", "sims_1b_extension_trt_mis_xsi0_bp.jpg"), plot = bp.xsi0)
ggsave(here::here("results", "sims_1b_extension_trt_mis_xsi1_bp.jpg"), plot = bp.xsi1)


#One Stage: ego-alter treatment-free model components mis-specified

study1b_extension_trt_fr_mis <- readRDS(here::here("results", "study1b_extension_trt_fr_mis.rds"))

names(study1b_extension_trt_fr_mis) <- c("AOpts_study1b_ext_trt_fr_mis", "YOpts_study1b_ext_trt_fr_mis", "Regret_study1b_ext_trt_fr_mis", "Psi_Dw_study1b_ext_trt_fr_mis", "Psi_M1_study1b_ext_trt_fr_mis", "Psi_M2_study1b_ext_trt_fr_mis", "Psi_M3_study1b_ext_trt_fr_mis", "Psi_M4_study1b_ext_trt_fr_mis", "Psi_M5_study1b_ext_trt_fr_mis", "Xsi_Dw_study1b_ext_trt_fr_mis", "Xsi_M1_study1b_ext_trt_fr_mis", "Xsi_M2_study1b_ext_trt_fr_mis", "Xsi_M3_study1b_ext_trt_fr_mis", "Xsi_M4_study1b_ext_trt_fr_mis", "Xsi_M5_study1b_ext_trt_fr_mis", "Mean_W_alters_study1b_ext_trt_fr_mis", "Mean_B_alters_study1b_ext_trt_fr_mis")

AOpts.df <- study1b_extension_trt_fr_mis$AOpts_study1b_ext_trt_fr_mis
YOpts.df <- study1b_extension_trt_fr_mis$YOpts_study1b_ext_trt_fr_mis
Regret.df <- study1b_extension_trt_fr_mis$Regret_study1b_ext_trt_fr_mis

Psi.Dw <- study1b_extension_trt_fr_mis$Psi_Dw_study1b_ext_trt_fr_mis
Psi.M1 <- study1b_extension_trt_fr_mis$Psi_M1_study1b_ext_trt_fr_mis
Psi.M2 <- study1b_extension_trt_fr_mis$Psi_M2_study1b_ext_trt_fr_mis
Psi.M3 <- study1b_extension_trt_fr_mis$Psi_M3_study1b_ext_trt_fr_mis
Psi.M4 <- study1b_extension_trt_fr_mis$Psi_M4_study1b_ext_trt_fr_mis
Psi.M5 <- study1b_extension_trt_fr_mis$Psi_M5_study1b_ext_trt_fr_mis

Xsi.Dw <- study1b_extension_trt_fr_mis$Xsi_Dw_study1b_ext_trt_fr_mis
Xsi.M1 <- study1b_extension_trt_fr_mis$Xsi_M1_study1b_ext_trt_fr_mis
Xsi.M2 <- study1b_extension_trt_fr_mis$Xsi_M2_study1b_ext_trt_fr_mis
Xsi.M3 <- study1b_extension_trt_fr_mis$Xsi_M3_study1b_ext_trt_fr_mis
Xsi.M4 <- study1b_extension_trt_fr_mis$Xsi_M4_study1b_ext_trt_fr_mis
Xsi.M5 <- study1b_extension_trt_fr_mis$Xsi_M5_study1b_ext_trt_fr_mis

#Provides numerical, tabular summary of optimal treatment rates and regrets across methods
plot_trt_rates <- na.omit(AOpts.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(TrtRate),
         mean = mean(TrtRate),
         sd = sd(TrtRate),
         empse = sqrt((sum((TrtRate-mean(TrtRate))^2))/(n-1)),
         .groups="keep")

plot_trt_rates

plot_outcomes <- na.omit(YOpts.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(Outcome),
         mean = mean(Outcome),
         sd = sd(Outcome),
         empse = sqrt((sum((Outcome-mean(Outcome))^2))/(n-1)),
         .groups="keep")

plot_outcomes

plot_regrets <- na.omit(Regret.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(Regret),
         mean = mean(Regret),
         sd = sd(Regret),
         empse = sqrt((sum((Regret-mean(Regret))^2))/(n-1)),
         .groups="keep")

plot_regrets

write.csv(plot_trt_rates, here::here("results", "sims_1b_extension_trt_fr_mis_trtRate_tbl.csv"))

write.csv(plot_outcomes, here::here("results", "sims_1b_extension_trt_fr_mis_outcome_tbl.csv"))

write.csv(plot_regrets, here::here("results", "sims_1b_extension_trt_fr_mis_regret_tbl.csv"))

AOpts.df$phys.grps <- as.factor(AOpts.df$phys.grps)
YOpts.df$phys.grps <- as.factor(YOpts.df$phys.grps)

Psi.Dw$phys.grps <- as.factor(Psi.Dw$phys.grps)
Psi.Dw$Method <- rep("ME1", nrow(Psi.Dw))

Psi.M1$phys.grps <- as.factor(Psi.M1$phys.grps)
Psi.M1$Method <- rep("ME2", nrow(Psi.M1))

Psi.M2$phys.grps <- as.factor(Psi.M2$phys.grps)
Psi.M2$Method <- rep("ME3", nrow(Psi.M2))

Psi.M3$phys.grps <- as.factor(Psi.M3$phys.grps)
Psi.M3$Method <- rep("ME4", nrow(Psi.M3))

Psi.M4$phys.grps <- as.factor(Psi.M4$phys.grps)
Psi.M4$Method <- rep("ME5", nrow(Psi.M4))

Psi.M5$phys.grps <- as.factor(Psi.M5$phys.grps)
Psi.M5$Method <- rep("ME6", nrow(Psi.M5))

Xsi.Dw$phys.grps <- as.factor(Xsi.Dw$phys.grps)
Xsi.Dw$Method <- rep("ME1", nrow(Xsi.Dw))

Xsi.M1$phys.grps <- as.factor(Xsi.M1$phys.grps)
Xsi.M1$Method <- rep("ME2", nrow(Xsi.M1))

Xsi.M2$phys.grps <- as.factor(Xsi.M2$phys.grps)
Xsi.M2$Method <- rep("ME3", nrow(Xsi.M2))

Xsi.M3$phys.grps <- as.factor(Xsi.M3$phys.grps)
Xsi.M3$Method <- rep("ME4", nrow(Xsi.M3))

Xsi.M4$phys.grps <- as.factor(Xsi.M4$phys.grps)
Xsi.M4$Method <- rep("ME5", nrow(Xsi.M4))

Xsi.M5$phys.grps <- as.factor(Xsi.M5$phys.grps)
Xsi.M5$Method <- rep("ME6", nrow(Xsi.M5))


Psi.0.df <- rbind(Psi.Dw[Psi.Dw$PsiParameter=="psi.0",], Psi.M1[Psi.M1$PsiParameter=="psi.0",], Psi.M2[Psi.M2$PsiParameter=="psi.0",], Psi.M3[Psi.M3$PsiParameter=="psi.0",], Psi.M4[Psi.M4$PsiParameter=="psi.0",], Psi.M5[Psi.M5$PsiParameter=="psi.0",])
Psi.1.df <- rbind(Psi.Dw[Psi.Dw$PsiParameter=="psi.1",], Psi.M1[Psi.M1$PsiParameter=="psi.1",], Psi.M2[Psi.M2$PsiParameter=="psi.1",], Psi.M3[Psi.M3$PsiParameter=="psi.1",], Psi.M4[Psi.M4$PsiParameter=="psi.1",], Psi.M5[Psi.M5$PsiParameter=="psi.1",])
Psi.2.df <- rbind(Psi.Dw[Psi.Dw$PsiParameter=="psi.2",], Psi.M1[Psi.M1$PsiParameter=="psi.2",], Psi.M2[Psi.M2$PsiParameter=="psi.2",], Psi.M3[Psi.M3$PsiParameter=="psi.2",], Psi.M4[Psi.M4$PsiParameter=="psi.2",], Psi.M5[Psi.M5$PsiParameter=="psi.2",])
Psi.3.df <- rbind(Psi.Dw[Psi.Dw$PsiParameter=="psi.3",], Psi.M1[Psi.M1$PsiParameter=="psi.3",], Psi.M2[Psi.M2$PsiParameter=="psi.3",], Psi.M3[Psi.M3$PsiParameter=="psi.3",], Psi.M4[Psi.M4$PsiParameter=="psi.3",], Psi.M5[Psi.M5$PsiParameter=="psi.3",])

Xsi.0.df <- rbind(Xsi.Dw[Xsi.Dw$XsiParameter=="xsi.0",], Xsi.M1[Xsi.M1$XsiParameter=="xsi.0",], Xsi.M2[Xsi.M2$XsiParameter=="xsi.0",], Xsi.M3[Xsi.M3$XsiParameter=="xsi.0",], Xsi.M4[Xsi.M4$XsiParameter=="xsi.0",], Xsi.M5[Xsi.M5$XsiParameter=="xsi.0",])
Xsi.1.df <- rbind(Xsi.Dw[Xsi.Dw$XsiParameter=="xsi.1",], Xsi.M1[Xsi.M1$XsiParameter=="xsi.1",], Xsi.M2[Xsi.M2$XsiParameter=="xsi.1",], Xsi.M3[Xsi.M3$XsiParameter=="xsi.1",], Xsi.M4[Xsi.M4$XsiParameter=="xsi.1",], Xsi.M5[Xsi.M5$XsiParameter=="xsi.1",])

bp.trt <- ggplot(AOpts.df, aes(x=Method, y=TrtRate, fill=Method)) + geom_boxplot() + labs(y="Optimal treatment rate") + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + theme(text = element_text(size = 14))
bp.out <- ggplot(YOpts.df, aes(x=Method, y=Outcome, fill=Method)) + geom_boxplot() + labs(y="Optimal outcome") + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + theme(text = element_text(size = 14))
bp.psi0 <- ggplot(Psi.0.df, aes(x=Method, y=Estimate, fill=Method)) + geom_boxplot() + labs(y=TeX(r"($\psi_{0}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi1 <- ggplot(Psi.1.df, aes(x=Method, y=Estimate, fill=Method)) + geom_boxplot() + labs(y=TeX(r"($\psi_{1}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi2 <- ggplot(Psi.2.df, aes(x=Method, y=Estimate, fill=Method)) + geom_boxplot() + labs(y=TeX(r"($\psi_{2}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=-1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi3 <- ggplot(Psi.3.df, aes(x=Method, y=Estimate, fill=Method)) + geom_boxplot() + labs(y=TeX(r"($\psi_{3}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=-1.5, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.xsi0 <- ggplot(Xsi.0.df, aes(x=Method, y=Estimate, fill=Method)) + geom_boxplot() + labs(y=TeX(r"($\xi_{0}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.xsi1 <- ggplot(Xsi.1.df, aes(x=Method, y=Estimate, fill=Method)) + geom_boxplot() + labs(y=TeX(r"($\xi_{1}$ Estimate)")) + scale_fill_brewer(palette = "Greys") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))

ggsave(here::here("results", "sims_1b_extension_trt_fr_mis_trtRate_bp.pdf"), plot = bp.trt)
ggsave(here::here("results", "sims_1b_extension_trt_fr_mis_outcome_bp.pdf"), plot = bp.out)
ggsave(here::here("results", "sims_1b_extension_trt_fr_mis_psi0_bp.pdf"), plot = bp.psi0)
ggsave(here::here("results", "sims_1b_extension_trt_fr_mis_psi1_bp.pdf"), plot = bp.psi1)
ggsave(here::here("results", "sims_1b_extension_trt_fr_mis_psi2_bp.pdf"), plot = bp.psi2)
ggsave(here::here("results", "sims_1b_extension_trt_fr_mis_psi3_bp.pdf"), plot = bp.psi3)
ggsave(here::here("results", "sims_1b_extension_trt_fr_mis_xsi0_bp.pdf"), plot = bp.xsi0)
ggsave(here::here("results", "sims_1b_extension_trt_fr_mis_xsi1_bp.pdf"), plot = bp.xsi1)

ggsave(here::here("results", "sims_1b_extension_trt_fr_mis_trtRate_bp.jpg"), plot = bp.trt)
ggsave(here::here("results", "sims_1b_extension_trt_fr_mis_outcome_bp.jpg"), plot = bp.out)
ggsave(here::here("results", "sims_1b_extension_trt_fr_mis_psi0_bp.jpg"), plot = bp.psi0)
ggsave(here::here("results", "sims_1b_extension_trt_fr_mis_psi1_bp.jpg"), plot = bp.psi1)
ggsave(here::here("results", "sims_1b_extension_trt_fr_mis_psi2_bp.jpg"), plot = bp.psi2)
ggsave(here::here("results", "sims_1b_extension_trt_fr_mis_psi3_bp.jpg"), plot = bp.psi3)
ggsave(here::here("results", "sims_1b_extension_trt_fr_mis_xsi0_bp.jpg"), plot = bp.xsi0)
ggsave(here::here("results", "sims_1b_extension_trt_fr_mis_xsi1_bp.jpg"), plot = bp.xsi1)



#***************************************************************************************************

study1b_1to5 <- readRDS(here::here("results", "study1b_1to5.rds"))

names(study1b_1to5) <- c("AOpts_study1b_1to5", "YOpts_study1b_1to5", "Regret_study1b_1to5", "Psi_Dw_study1b_1to5", "Psi_M1_study1b_1to5", "Psi_M2_study1b_1to5", "Psi_M3_study1b_1to5", "Psi_M4_study1b_1to5", "Psi_M5_study1b_1to5", "Psi_M6_study1b_1to5", "Psi_M7_study1b_1to5", "Psi_M8_study1b_1to5", "Psi_M9_study1b_1to5", "Xsi_Dw_study1b_1to5", "Xsi_M1_study1b_1to5", "Xsi_M2_study1b_1to5", "Xsi_M3_study1b_1to5", "Xsi_M4_study1b_1to5", "Xsi_M5_study1b_1to5", "Xsi_M6_study1b_1to5", "Xsi_M7_study1b_1to5", "Xsi_M8_study1b_1to5", "Xsi_M9_study1b_1to5", "Mean_W_alters_study1b_1to5", "Mean_B_alters_study1b_1to5")

#One Stage: ego-alter treatment-free model components mis-specified, treatment model correct, scenarios 1-5

AOpts.df <- study1b_1to5$AOpts_study1b_1to5
YOpts.df <- study1b_1to5$YOpts_study1b_1to5
Regret.df <- study1b_1to5$Regret_study1b_1to5

Psi.Dw <- study1b_1to5$Psi_Dw_study1b_1to5
Psi.M1 <- study1b_1to5$Psi_M1_study1b_1to5
Psi.M2 <- study1b_1to5$Psi_M2_study1b_1to5
Psi.M3 <- study1b_1to5$Psi_M3_study1b_1to5
Psi.M4 <- study1b_1to5$Psi_M4_study1b_1to5
Psi.M5 <- study1b_1to5$Psi_M5_study1b_1to5
Psi.M6 <- study1b_1to5$Psi_M6_study1b_1to5
Psi.M7 <- study1b_1to5$Psi_M7_study1b_1to5
Psi.M8 <- study1b_1to5$Psi_M8_study1b_1to5
Psi.M9 <- study1b_1to5$Psi_M9_study1b_1to5

Xsi.Dw <- study1b_1to5$Xsi_Dw_study1b_1to5
Xsi.M1 <- study1b_1to5$Xsi_M1_study1b_1to5
Xsi.M2 <- study1b_1to5$Xsi_M2_study1b_1to5
Xsi.M3 <- study1b_1to5$Xsi_M3_study1b_1to5
Xsi.M4 <- study1b_1to5$Xsi_M4_study1b_1to5
Xsi.M5 <- study1b_1to5$Xsi_M5_study1b_1to5
Xsi.M6 <- study1b_1to5$Xsi_M6_study1b_1to5
Xsi.M7 <- study1b_1to5$Xsi_M7_study1b_1to5
Xsi.M8 <- study1b_1to5$Xsi_M8_study1b_1to5
Xsi.M9 <- study1b_1to5$Xsi_M9_study1b_1to5

#Provides numerical, tabular summary of optimal treatment rates and regrets across methods
plot_trt_rates <- na.omit(AOpts.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(TrtRate),
         mean = mean(TrtRate),
         sd = sd(TrtRate),
         empse = sqrt((sum((TrtRate-mean(TrtRate))^2))/(n-1)),
         .groups="keep")

plot_trt_rates

plot_outcomes <- na.omit(YOpts.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(Outcome),
         mean = mean(Outcome),
         sd = sd(Outcome),
         empse = sqrt((sum((Outcome-mean(Outcome))^2))/(n-1)),
         .groups="keep")

plot_outcomes

plot_regrets <- na.omit(Regret.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(Regret),
         mean = mean(Regret),
         sd = sd(Regret),
         empse = sqrt((sum((Regret-mean(Regret))^2))/(n-1)),
         .groups="keep")

plot_regrets

write.csv(plot_trt_rates, here::here("results", "sims_1b_1to5_trtRate_tbl.csv"))

write.csv(plot_outcomes, here::here("results", "sims_1b_1to5_outcome_tbl.csv"))

write.csv(plot_regrets, here::here("results", "sims_1b_1to5_regret_tbl.csv"))

AOpts.df$phys.grps <- as.factor(AOpts.df$phys.grps)
YOpts.df$phys.grps <- as.factor(YOpts.df$phys.grps)

Psi.Dw$phys.grps <- as.factor(Psi.Dw$phys.grps)
Psi.Dw$Method <- rep("ME1", nrow(Psi.Dw))

Psi.M1$phys.grps <- as.factor(Psi.M1$phys.grps)
Psi.M1$Method <- rep("ME2", nrow(Psi.M1))

Psi.M2$phys.grps <- as.factor(Psi.M2$phys.grps)
Psi.M2$Method <- rep("ME3", nrow(Psi.M2))

Psi.M3$phys.grps <- as.factor(Psi.M3$phys.grps)
Psi.M3$Method <- rep("ME4", nrow(Psi.M3))

Psi.M4$phys.grps <- as.factor(Psi.M4$phys.grps)
Psi.M4$Method <- rep("ME5", nrow(Psi.M4))

Psi.M5$phys.grps <- as.factor(Psi.M5$phys.grps)
Psi.M5$Method <- rep("ME6", nrow(Psi.M5))

Psi.M6$phys.grps <- as.factor(Psi.M6$phys.grps)
Psi.M6$Method <- rep("ME7", nrow(Psi.M6))

Psi.M7$phys.grps <- as.factor(Psi.M7$phys.grps)
Psi.M7$Method <- rep("ME8", nrow(Psi.M7))

Psi.M8$phys.grps <- as.factor(Psi.M8$phys.grps)
Psi.M8$Method <- rep("ME9", nrow(Psi.M8))

Psi.M9$phys.grps <- as.factor(Psi.M9$phys.grps)
Psi.M9$Method <- rep("ME10", nrow(Psi.M9))

Xsi.Dw$phys.grps <- as.factor(Xsi.Dw$phys.grps)
Xsi.Dw$Method <- rep("ME1", nrow(Xsi.Dw))

Xsi.M1$phys.grps <- as.factor(Xsi.M1$phys.grps)
Xsi.M1$Method <- rep("ME2", nrow(Xsi.M1))

Xsi.M2$phys.grps <- as.factor(Xsi.M2$phys.grps)
Xsi.M2$Method <- rep("ME3", nrow(Xsi.M2))

Xsi.M3$phys.grps <- as.factor(Xsi.M3$phys.grps)
Xsi.M3$Method <- rep("ME4", nrow(Xsi.M3))

Xsi.M4$phys.grps <- as.factor(Xsi.M4$phys.grps)
Xsi.M4$Method <- rep("ME5", nrow(Xsi.M4))

Xsi.M5$phys.grps <- as.factor(Xsi.M5$phys.grps)
Xsi.M5$Method <- rep("ME6", nrow(Xsi.M5))

Xsi.M6$phys.grps <- as.factor(Xsi.M6$phys.grps)
Xsi.M6$Method <- rep("ME7", nrow(Xsi.M6))

Xsi.M7$phys.grps <- as.factor(Xsi.M7$phys.grps)
Xsi.M7$Method <- rep("ME8", nrow(Xsi.M7))

Xsi.M8$phys.grps <- as.factor(Xsi.M8$phys.grps)
Xsi.M8$Method <- rep("ME9", nrow(Xsi.M8))

Xsi.M9$phys.grps <- as.factor(Xsi.M9$phys.grps)
Xsi.M9$Method <- rep("ME10", nrow(Xsi.M9))

Psi.0.df <- rbind(Psi.Dw[Psi.Dw$PsiParameter=="psi.0",], Psi.M1[Psi.M1$PsiParameter=="psi.0",], Psi.M2[Psi.M2$PsiParameter=="psi.0",], Psi.M3[Psi.M3$PsiParameter=="psi.0",], Psi.M4[Psi.M4$PsiParameter=="psi.0",], Psi.M5[Psi.M5$PsiParameter=="psi.0",], Psi.M6[Psi.M6$PsiParameter=="psi.0",], Psi.M7[Psi.M7$PsiParameter=="psi.0",], Psi.M8[Psi.M8$PsiParameter=="psi.0",], Psi.M9[Psi.M9$PsiParameter=="psi.0",])
Psi.1.df <- rbind(Psi.Dw[Psi.Dw$PsiParameter=="psi.1",], Psi.M1[Psi.M1$PsiParameter=="psi.1",], Psi.M2[Psi.M2$PsiParameter=="psi.1",], Psi.M3[Psi.M3$PsiParameter=="psi.1",], Psi.M4[Psi.M4$PsiParameter=="psi.1",], Psi.M5[Psi.M5$PsiParameter=="psi.1",], Psi.M6[Psi.M6$PsiParameter=="psi.1",], Psi.M7[Psi.M7$PsiParameter=="psi.1",], Psi.M8[Psi.M8$PsiParameter=="psi.1",], Psi.M9[Psi.M9$PsiParameter=="psi.1",])
Psi.2.df <- rbind(Psi.Dw[Psi.Dw$PsiParameter=="psi.2",], Psi.M1[Psi.M1$PsiParameter=="psi.2",], Psi.M2[Psi.M2$PsiParameter=="psi.2",], Psi.M3[Psi.M3$PsiParameter=="psi.2",], Psi.M4[Psi.M4$PsiParameter=="psi.2",], Psi.M5[Psi.M5$PsiParameter=="psi.2",], Psi.M6[Psi.M6$PsiParameter=="psi.2",], Psi.M7[Psi.M7$PsiParameter=="psi.2",], Psi.M8[Psi.M8$PsiParameter=="psi.2",], Psi.M9[Psi.M9$PsiParameter=="psi.2",])
Psi.3.df <- rbind(Psi.Dw[Psi.Dw$PsiParameter=="psi.3",], Psi.M1[Psi.M1$PsiParameter=="psi.3",], Psi.M2[Psi.M2$PsiParameter=="psi.3",], Psi.M3[Psi.M3$PsiParameter=="psi.3",], Psi.M4[Psi.M4$PsiParameter=="psi.3",], Psi.M5[Psi.M5$PsiParameter=="psi.3",], Psi.M6[Psi.M6$PsiParameter=="psi.3",], Psi.M7[Psi.M7$PsiParameter=="psi.3",], Psi.M8[Psi.M8$PsiParameter=="psi.3",], Psi.M9[Psi.M9$PsiParameter=="psi.3",])

Xsi.0.df <- rbind(Xsi.Dw[Xsi.Dw$XsiParameter=="xsi.0",], Xsi.M1[Xsi.M1$XsiParameter=="xsi.0",], Xsi.M2[Xsi.M2$XsiParameter=="xsi.0",], Xsi.M3[Xsi.M3$XsiParameter=="xsi.0",], Xsi.M4[Xsi.M4$XsiParameter=="xsi.0",], Xsi.M5[Xsi.M5$XsiParameter=="xsi.0",], Xsi.M6[Xsi.M6$XsiParameter=="xsi.0",], Xsi.M7[Xsi.M7$XsiParameter=="xsi.0",], Xsi.M8[Xsi.M8$XsiParameter=="xsi.0",], Xsi.M9[Xsi.M9$XsiParameter=="xsi.0",])
Xsi.1.df <- rbind(Xsi.Dw[Xsi.Dw$XsiParameter=="xsi.1",], Xsi.M1[Xsi.M1$XsiParameter=="xsi.1",], Xsi.M2[Xsi.M2$XsiParameter=="xsi.1",], Xsi.M3[Xsi.M3$XsiParameter=="xsi.1",], Xsi.M4[Xsi.M4$XsiParameter=="xsi.1",], Xsi.M5[Xsi.M5$XsiParameter=="xsi.1",], Xsi.M6[Xsi.M6$XsiParameter=="xsi.1",], Xsi.M7[Xsi.M7$XsiParameter=="xsi.1",], Xsi.M8[Xsi.M8$XsiParameter=="xsi.1",], Xsi.M9[Xsi.M9$XsiParameter=="xsi.1",])

level_order_ME <- c("Observed", "ME1", "ME2", "ME3", "ME4", "ME5", "ME6", "ME7", "ME8", "ME9", "ME10")

bp.trt <- ggplot(AOpts.df, aes(x=factor(Method, level=level_order_ME), y=TrtRate, fill=Method)) + geom_boxplot(fill="grey") + labs(y="Optimal treatment rate") + labs(x="Scenario") + theme_classic() + theme(legend.position = "none") + theme(text = element_text(size = 14))
bp.out <- ggplot(YOpts.df, aes(x=factor(Method, level=level_order_ME), y=Outcome, fill=Method)) + geom_boxplot(fill="grey") + labs(y="Optimal outcome") + labs(x="Scenario") + theme_classic() + theme(legend.position = "none") + theme(text = element_text(size = 14))
bp.psi0 <- ggplot(Psi.0.df, aes(x=factor(Method, level=level_order_ME), y=Estimate, fill=Method)) + geom_boxplot(fill="grey") + labs(y=TeX(r"($\psi_{0}$ Estimate)")) + labs(x="Scenario") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi1 <- ggplot(Psi.1.df, aes(x=factor(Method, level=level_order_ME), y=Estimate, fill=Method)) + geom_boxplot(fill="grey") + labs(y=TeX(r"($\psi_{1}$ Estimate)")) + labs(x="Scenario") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi2 <- ggplot(Psi.2.df, aes(x=factor(Method, level=level_order_ME), y=Estimate, fill=Method)) + geom_boxplot(fill="grey") + labs(y=TeX(r"($\psi_{2}$ Estimate)")) + labs(x="Scenario") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=-1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi3 <- ggplot(Psi.3.df, aes(x=factor(Method, level=level_order_ME), y=Estimate, fill=Method)) + geom_boxplot(fill="grey") + labs(y=TeX(r"($\psi_{3}$ Estimate)")) + labs(x="Scenario") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=-1.5, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.xsi0 <- ggplot(Xsi.0.df, aes(x=factor(Method, level=level_order_ME), y=Estimate, fill=Method)) + geom_boxplot(fill="grey") + labs(y=TeX(r"($\xi_{0}$ Estimate)")) + labs(x="Scenario") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.xsi1 <- ggplot(Xsi.1.df, aes(x=factor(Method, level=level_order_ME), y=Estimate, fill=Method)) + geom_boxplot(fill="grey") + labs(y=TeX(r"($\xi_{1}$ Estimate)")) + labs(x="Scenario") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))

ggsave(here::here("results", "sims_1b_1to5_trtRate_bp.pdf"), plot = bp.trt)
ggsave(here::here("results", "sims_1b_1to5_outcome_bp.pdf"), plot = bp.out)
ggsave(here::here("results", "sims_1b_1to5_psi0_bp.pdf"), plot = bp.psi0)
ggsave(here::here("results", "sims_1b_1to5_psi1_bp.pdf"), plot = bp.psi1)
ggsave(here::here("results", "sims_1b_1to5_psi2_bp.pdf"), plot = bp.psi2)
ggsave(here::here("results", "sims_1b_1to5_psi3_bp.pdf"), plot = bp.psi3)
ggsave(here::here("results", "sims_1b_1to5_xsi0_bp.pdf"), plot = bp.xsi0)
ggsave(here::here("results", "sims_1b_1to5_xsi1_bp.pdf"), plot = bp.xsi1)

ggsave(here::here("results", "sims_1b_1to5_trtRate_bp.jpg"), plot = bp.trt)
ggsave(here::here("results", "sims_1b_1to5_outcome_bp.jpg"), plot = bp.out)
ggsave(here::here("results", "sims_1b_1to5_psi0_bp.jpg"), plot = bp.psi0)
ggsave(here::here("results", "sims_1b_1to5_psi1_bp.jpg"), plot = bp.psi1)
ggsave(here::here("results", "sims_1b_1to5_psi2_bp.jpg"), plot = bp.psi2)
ggsave(here::here("results", "sims_1b_1to5_psi3_bp.jpg"), plot = bp.psi3)
ggsave(here::here("results", "sims_1b_1to5_xsi0_bp.jpg"), plot = bp.xsi0)
ggsave(here::here("results", "sims_1b_1to5_xsi1_bp.jpg"), plot = bp.xsi1)


#***************************************************************************************************


#One Stage: ego-alter treatment-free model components mis-specified, treatment model incorrect, scenarios 6-10

study1b_6to10 <- readRDS(here::here("results", "study1b_6to10.rds"))

names(study1b_6to10) <- c("AOpts_study1b_6to10", "YOpts_study1b_6to10", "Regret_study1b_6to10", "Psi_Dw_study1b_6to10", "Psi_M1_study1b_6to10", "Psi_M2_study1b_6to10", "Psi_M3_study1b_6to10", "Psi_M4_study1b_6to10", "Psi_M5_study1b_6to10", "Psi_M6_study1b_6to10", "Psi_M7_study1b_6to10", "Psi_M8_study1b_6to10", "Psi_M9_study1b_6to10", "Xsi_Dw_study1b_6to10", "Xsi_M1_study1b_6to10", "Xsi_M2_study1b_6to10", "Xsi_M3_study1b_6to10", "Xsi_M4_study1b_6to10", "Xsi_M5_study1b_6to10", "Xsi_M6_study1b_6to10", "Xsi_M7_study1b_6to10", "Xsi_M8_study1b_6to10", "Xsi_M9_study1b_6to10", "Mean_W_alters_study1b_6to10", "Mean_B_alters_study1b_6to10")

AOpts.df <- study1b_6to10$AOpts_study1b_6to10
YOpts.df <- study1b_6to10$YOpts_study1b_6to10
Regret.df <- study1b_6to10$Regret_study1b_6to10

Psi.Dw <- study1b_6to10$Psi_Dw_study1b_6to10
Psi.M1 <- study1b_6to10$Psi_M1_study1b_6to10
Psi.M2 <- study1b_6to10$Psi_M2_study1b_6to10
Psi.M3 <- study1b_6to10$Psi_M3_study1b_6to10
Psi.M4 <- study1b_6to10$Psi_M4_study1b_6to10
Psi.M5 <- study1b_6to10$Psi_M5_study1b_6to10
Psi.M6 <- study1b_6to10$Psi_M6_study1b_6to10
Psi.M7 <- study1b_6to10$Psi_M7_study1b_6to10
Psi.M8 <- study1b_6to10$Psi_M8_study1b_6to10
Psi.M9 <- study1b_6to10$Psi_M9_study1b_6to10

Xsi.Dw <- study1b_6to10$Xsi_Dw_study1b_6to10
Xsi.M1 <- study1b_6to10$Xsi_M1_study1b_6to10
Xsi.M2 <- study1b_6to10$Xsi_M2_study1b_6to10
Xsi.M3 <- study1b_6to10$Xsi_M3_study1b_6to10
Xsi.M4 <- study1b_6to10$Xsi_M4_study1b_6to10
Xsi.M5 <- study1b_6to10$Xsi_M5_study1b_6to10
Xsi.M6 <- study1b_6to10$Xsi_M6_study1b_6to10
Xsi.M7 <- study1b_6to10$Xsi_M7_study1b_6to10
Xsi.M8 <- study1b_6to10$Xsi_M8_study1b_6to10
Xsi.M9 <- study1b_6to10$Xsi_M9_study1b_6to10

#Provides numerical, tabular summary of optimal treatment rates and regrets across methods
plot_trt_rates <- na.omit(AOpts.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(TrtRate),
         mean = mean(TrtRate),
         sd = sd(TrtRate),
         empse = sqrt((sum((TrtRate-mean(TrtRate))^2))/(n-1)),
         .groups="keep")

plot_trt_rates

plot_outcomes <- na.omit(YOpts.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(Outcome),
         mean = mean(Outcome),
         sd = sd(Outcome),
         empse = sqrt((sum((Outcome-mean(Outcome))^2))/(n-1)),
         .groups="keep")

plot_outcomes

plot_regrets <- na.omit(Regret.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(Regret),
         mean = mean(Regret),
         sd = sd(Regret),
         empse = sqrt((sum((Regret-mean(Regret))^2))/(n-1)),
         .groups="keep")

plot_regrets

write.csv(plot_trt_rates, here::here("results", "sims_1b_6to10_trtRate_tbl.csv"))

write.csv(plot_outcomes, here::here("results", "sims_1b_6to10_outcome_tbl.csv"))

write.csv(plot_regrets, here::here("results", "sims_1b_6to10_regret_tbl.csv"))

AOpts.df$phys.grps <- as.factor(AOpts.df$phys.grps)
YOpts.df$phys.grps <- as.factor(YOpts.df$phys.grps)

Psi.Dw$phys.grps <- as.factor(Psi.Dw$phys.grps)
Psi.Dw$Method <- rep("ME1", nrow(Psi.Dw))

Psi.M1$phys.grps <- as.factor(Psi.M1$phys.grps)
Psi.M1$Method <- rep("ME2", nrow(Psi.M1))

Psi.M2$phys.grps <- as.factor(Psi.M2$phys.grps)
Psi.M2$Method <- rep("ME3", nrow(Psi.M2))

Psi.M3$phys.grps <- as.factor(Psi.M3$phys.grps)
Psi.M3$Method <- rep("ME4", nrow(Psi.M3))

Psi.M4$phys.grps <- as.factor(Psi.M4$phys.grps)
Psi.M4$Method <- rep("ME5", nrow(Psi.M4))

Psi.M5$phys.grps <- as.factor(Psi.M5$phys.grps)
Psi.M5$Method <- rep("ME6", nrow(Psi.M5))

Psi.M6$phys.grps <- as.factor(Psi.M6$phys.grps)
Psi.M6$Method <- rep("ME7", nrow(Psi.M6))

Psi.M7$phys.grps <- as.factor(Psi.M7$phys.grps)
Psi.M7$Method <- rep("ME8", nrow(Psi.M7))

Psi.M8$phys.grps <- as.factor(Psi.M8$phys.grps)
Psi.M8$Method <- rep("ME9", nrow(Psi.M8))

Psi.M9$phys.grps <- as.factor(Psi.M9$phys.grps)
Psi.M9$Method <- rep("ME10", nrow(Psi.M9))

Xsi.Dw$phys.grps <- as.factor(Xsi.Dw$phys.grps)
Xsi.Dw$Method <- rep("ME1", nrow(Xsi.Dw))

Xsi.M1$phys.grps <- as.factor(Xsi.M1$phys.grps)
Xsi.M1$Method <- rep("ME2", nrow(Xsi.M1))

Xsi.M2$phys.grps <- as.factor(Xsi.M2$phys.grps)
Xsi.M2$Method <- rep("ME3", nrow(Xsi.M2))

Xsi.M3$phys.grps <- as.factor(Xsi.M3$phys.grps)
Xsi.M3$Method <- rep("ME4", nrow(Xsi.M3))

Xsi.M4$phys.grps <- as.factor(Xsi.M4$phys.grps)
Xsi.M4$Method <- rep("ME5", nrow(Xsi.M4))

Xsi.M5$phys.grps <- as.factor(Xsi.M5$phys.grps)
Xsi.M5$Method <- rep("ME6", nrow(Xsi.M5))

Xsi.M6$phys.grps <- as.factor(Xsi.M6$phys.grps)
Xsi.M6$Method <- rep("ME7", nrow(Xsi.M6))

Xsi.M7$phys.grps <- as.factor(Xsi.M7$phys.grps)
Xsi.M7$Method <- rep("ME8", nrow(Xsi.M7))

Xsi.M8$phys.grps <- as.factor(Xsi.M8$phys.grps)
Xsi.M8$Method <- rep("ME9", nrow(Xsi.M8))

Xsi.M9$phys.grps <- as.factor(Xsi.M9$phys.grps)
Xsi.M9$Method <- rep("ME10", nrow(Xsi.M9))

Psi.0.df <- rbind(Psi.Dw[Psi.Dw$PsiParameter=="psi.0",], Psi.M1[Psi.M1$PsiParameter=="psi.0",], Psi.M2[Psi.M2$PsiParameter=="psi.0",], Psi.M3[Psi.M3$PsiParameter=="psi.0",], Psi.M4[Psi.M4$PsiParameter=="psi.0",], Psi.M5[Psi.M5$PsiParameter=="psi.0",], Psi.M6[Psi.M6$PsiParameter=="psi.0",], Psi.M7[Psi.M7$PsiParameter=="psi.0",], Psi.M8[Psi.M8$PsiParameter=="psi.0",], Psi.M9[Psi.M9$PsiParameter=="psi.0",])
Psi.1.df <- rbind(Psi.Dw[Psi.Dw$PsiParameter=="psi.1",], Psi.M1[Psi.M1$PsiParameter=="psi.1",], Psi.M2[Psi.M2$PsiParameter=="psi.1",], Psi.M3[Psi.M3$PsiParameter=="psi.1",], Psi.M4[Psi.M4$PsiParameter=="psi.1",], Psi.M5[Psi.M5$PsiParameter=="psi.1",], Psi.M6[Psi.M6$PsiParameter=="psi.1",], Psi.M7[Psi.M7$PsiParameter=="psi.1",], Psi.M8[Psi.M8$PsiParameter=="psi.1",], Psi.M9[Psi.M9$PsiParameter=="psi.1",])
Psi.2.df <- rbind(Psi.Dw[Psi.Dw$PsiParameter=="psi.2",], Psi.M1[Psi.M1$PsiParameter=="psi.2",], Psi.M2[Psi.M2$PsiParameter=="psi.2",], Psi.M3[Psi.M3$PsiParameter=="psi.2",], Psi.M4[Psi.M4$PsiParameter=="psi.2",], Psi.M5[Psi.M5$PsiParameter=="psi.2",], Psi.M6[Psi.M6$PsiParameter=="psi.2",], Psi.M7[Psi.M7$PsiParameter=="psi.2",], Psi.M8[Psi.M8$PsiParameter=="psi.2",], Psi.M9[Psi.M9$PsiParameter=="psi.2",])
Psi.3.df <- rbind(Psi.Dw[Psi.Dw$PsiParameter=="psi.3",], Psi.M1[Psi.M1$PsiParameter=="psi.3",], Psi.M2[Psi.M2$PsiParameter=="psi.3",], Psi.M3[Psi.M3$PsiParameter=="psi.3",], Psi.M4[Psi.M4$PsiParameter=="psi.3",], Psi.M5[Psi.M5$PsiParameter=="psi.3",], Psi.M6[Psi.M6$PsiParameter=="psi.3",], Psi.M7[Psi.M7$PsiParameter=="psi.3",], Psi.M8[Psi.M8$PsiParameter=="psi.3",], Psi.M9[Psi.M9$PsiParameter=="psi.3",])

Xsi.0.df <- rbind(Xsi.Dw[Xsi.Dw$XsiParameter=="xsi.0",], Xsi.M1[Xsi.M1$XsiParameter=="xsi.0",], Xsi.M2[Xsi.M2$XsiParameter=="xsi.0",], Xsi.M3[Xsi.M3$XsiParameter=="xsi.0",], Xsi.M4[Xsi.M4$XsiParameter=="xsi.0",], Xsi.M5[Xsi.M5$XsiParameter=="xsi.0",], Xsi.M6[Xsi.M6$XsiParameter=="xsi.0",], Xsi.M7[Xsi.M7$XsiParameter=="xsi.0",], Xsi.M8[Xsi.M8$XsiParameter=="xsi.0",], Xsi.M9[Xsi.M9$XsiParameter=="xsi.0",])
Xsi.1.df <- rbind(Xsi.Dw[Xsi.Dw$XsiParameter=="xsi.1",], Xsi.M1[Xsi.M1$XsiParameter=="xsi.1",], Xsi.M2[Xsi.M2$XsiParameter=="xsi.1",], Xsi.M3[Xsi.M3$XsiParameter=="xsi.1",], Xsi.M4[Xsi.M4$XsiParameter=="xsi.1",], Xsi.M5[Xsi.M5$XsiParameter=="xsi.1",], Xsi.M6[Xsi.M6$XsiParameter=="xsi.1",], Xsi.M7[Xsi.M7$XsiParameter=="xsi.1",], Xsi.M8[Xsi.M8$XsiParameter=="xsi.1",], Xsi.M9[Xsi.M9$XsiParameter=="xsi.1",])

level_order_ME <- c("Observed", "ME1", "ME2", "ME3", "ME4", "ME5", "ME6", "ME7", "ME8", "ME9", "ME10")

bp.trt <- ggplot(AOpts.df, aes(x=factor(Method, level=level_order_ME), y=TrtRate, fill=Method)) + geom_boxplot(fill="grey") + labs(y="Optimal treatment rate") + labs(x="Scenario") + theme_classic() + theme(legend.position = "none") + theme(text = element_text(size = 14))
bp.out <- ggplot(YOpts.df, aes(x=factor(Method, level=level_order_ME), y=Outcome, fill=Method)) + geom_boxplot(fill="grey") + labs(y="Optimal outcome") + labs(x="Scenario") + theme_classic() + theme(legend.position = "none") + theme(text = element_text(size = 14))
bp.psi0 <- ggplot(Psi.0.df, aes(x=factor(Method, level=level_order_ME), y=Estimate, fill=Method)) + geom_boxplot(fill="grey") + labs(y=TeX(r"($\psi_{0}$ Estimate)")) + labs(x="Scenario") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi1 <- ggplot(Psi.1.df, aes(x=factor(Method, level=level_order_ME), y=Estimate, fill=Method)) + geom_boxplot(fill="grey") + labs(y=TeX(r"($\psi_{1}$ Estimate)")) + labs(x="Scenario") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi2 <- ggplot(Psi.2.df, aes(x=factor(Method, level=level_order_ME), y=Estimate, fill=Method)) + geom_boxplot(fill="grey") + labs(y=TeX(r"($\psi_{2}$ Estimate)")) + labs(x="Scenario") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=-1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.psi3 <- ggplot(Psi.3.df, aes(x=factor(Method, level=level_order_ME), y=Estimate, fill=Method)) + geom_boxplot(fill="grey") + labs(y=TeX(r"($\psi_{3}$ Estimate)")) + labs(x="Scenario") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=-1.5, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.xsi0 <- ggplot(Xsi.0.df, aes(x=factor(Method, level=level_order_ME), y=Estimate, fill=Method)) + geom_boxplot(fill="grey") + labs(y=TeX(r"($\xi_{0}$ Estimate)")) + labs(x="Scenario") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))
bp.xsi1 <- ggplot(Xsi.1.df, aes(x=factor(Method, level=level_order_ME), y=Estimate, fill=Method)) + geom_boxplot(fill="grey") + labs(y=TeX(r"($\xi_{1}$ Estimate)")) + labs(x="Scenario") + theme_classic() + theme(legend.position = "none") + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme(text = element_text(size = 14))

ggsave(here::here("results", "sims_1b_6to10_trtRate_bp.pdf"), plot = bp.trt)
ggsave(here::here("results", "sims_1b_6to10_outcome_bp.pdf"), plot = bp.out)
ggsave(here::here("results", "sims_1b_6to10_psi0_bp.pdf"), plot = bp.psi0)
ggsave(here::here("results", "sims_1b_6to10_psi1_bp.pdf"), plot = bp.psi1)
ggsave(here::here("results", "sims_1b_6to10_psi2_bp.pdf"), plot = bp.psi2)
ggsave(here::here("results", "sims_1b_6to10_psi3_bp.pdf"), plot = bp.psi3)
ggsave(here::here("results", "sims_1b_6to10_xsi0_bp.pdf"), plot = bp.xsi0)
ggsave(here::here("results", "sims_1b_6to10_xsi1_bp.pdf"), plot = bp.xsi1)

ggsave(here::here("results", "sims_1b_6to10_trtRate_bp.jpg"), plot = bp.trt)
ggsave(here::here("results", "sims_1b_6to10_outcome_bp.jpg"), plot = bp.out)
ggsave(here::here("results", "sims_1b_6to10_psi0_bp.jpg"), plot = bp.psi0)
ggsave(here::here("results", "sims_1b_6to10_psi1_bp.jpg"), plot = bp.psi1)
ggsave(here::here("results", "sims_1b_6to10_psi2_bp.jpg"), plot = bp.psi2)
ggsave(here::here("results", "sims_1b_6to10_psi3_bp.jpg"), plot = bp.psi3)
ggsave(here::here("results", "sims_1b_6to10_xsi0_bp.jpg"), plot = bp.xsi0)
ggsave(here::here("results", "sims_1b_6to10_xsi1_bp.jpg"), plot = bp.xsi1)


