#**********MAIN FILE: RUN DWOLS FOR WITHIN- AND BETWEEN-GROUP INTERFERENCE: ONE-STAGE SETTING**********

source(here::here("scripts", "sims_one_stage_trt_and_trt_fr_mis.R")) 

run <- 1000
n <- 50
sam.size <- 50
a <- 2
b <- 4
c <- 10
tau3 <- 1
tau2 <- 1
sigma2 <- 1
alpha.0 <- -1
alpha.1 <- 1
alpha.2 <- 0.5
psi.0 <- 1
psi.1 <- 1
psi.2 <- -1
psi.3 <- -1.5

bp_list <- begin.Sims(run, n=50, sam.size = 50, a=2, b=4, c=10, tau3=1, tau2=1, sigma2=1, alpha.0=-1, alpha.1=1, alpha.2=0.5, psi.0=1, psi.1=1, psi.2=-1, psi.3=-1.5)

J.list <- seq(a, b, 1)

phys.grps <- c()

phys.grps.regret <- c()

for(j in J.list){
  phys.val <- as.character(j)
  #Times five because we have what is observed plus four methods to compare
  phys.grps <- c(phys.grps, rep(phys.val, run*5))

  #For regrets, times four because we are comparing four methods
  phys.grps.regret <- c(phys.grps.regret, rep(phys.val, run*4))
}

AOpts.df <- bp_list$AOpts
AOpts.df$phys.grps <- phys.grps

YOpts.df <- bp_list$YOpts
YOpts.df$phys.grps <- phys.grps

YAgree.df <- bp_list$YAgree
YAgree.df$phys.grps <- phys.grps

Regret.df <- bp_list$RegretMethods
Regret.df$phys.grps <- phys.grps.regret

phys.dw1 <- c()

for(j in J.list){
  phys.val <- as.character(j)
  #Times 4 because we have four blip parameters
  phys.dw1 <- c(phys.dw1, rep(phys.val, run*4))
}

phys.dw2 <- c()

for(j in J.list){
  phys.val <- as.character(j)
  #Times 4 because we have four blip parameters
  phys.dw2 <- c(phys.dw2, rep(phys.val, run*4))
}

phys.dw3 <- c()

for(j in J.list){
  phys.val <- as.character(j)
  phys.dw3 <- c(phys.dw3, rep(phys.val, run*4))
}

phys.dw.xsi <- c()

for(j in J.list){
    phys.val <- as.character(j)
    phys.dw.xsi <- c(phys.dw.xsi, rep(phys.val, run*2))
}

Psi.Dw <- bp_list$Psi.Dw
Psi.Dw$phys.grps <- phys.dw1
Psi.Dw$TruePsi <- rep(rep(c(1, 1, -1, -1.5), each=run), length(J.list))

Psi.M1 <- bp_list$Psi.M1
Psi.M1$phys.grps <- phys.dw1
#"Each" is based off simulation runs
Psi.M1$TruePsi <- rep(rep(c(1, 1, -1, -1.5), each=run), length(J.list))

Psi.M2 <- bp_list$Psi.M2
Psi.M2$phys.grps <- phys.dw2
#"Each" is based off simulation runs
Psi.M2$TruePsi <- rep(rep(c(1, 1, -1, -1.5), each=run), length(J.list))

Psi.M3 <- bp_list$Psi.M3
Psi.M3$phys.grps <- phys.dw3
Psi.M3$TruePsi <- rep(rep(c(1, 1, -1, -1.5), each=run), length(J.list))

Xsi.Dw <- bp_list$Xsi.Dw
Xsi.Dw$phys.grps <- phys.dw.xsi
#"Each" is based off simulation runs
Xsi.Dw$TrueXsi <- rep(rep(c(1, 1), each=run), length(J.list))

Xsi.M1 <- bp_list$Xsi.M1
Xsi.M1$phys.grps <- phys.dw.xsi
#"Each" is based off simulation runs
Xsi.M1$TrueXsi <- rep(rep(c(1, 1), each=run), length(J.list))

Xsi.M2 <- bp_list$Xsi.M2
Xsi.M2$phys.grps <- phys.dw.xsi
#"Each" is based off simulation runs
Xsi.M2$TrueXsi <- rep(rep(c(1, 1), each=run), length(J.list))

Xsi.M3 <- bp_list$Xsi.M3
Xsi.M3$phys.grps <- phys.dw.xsi
Xsi.M3$TrueXsi <- rep(rep(c(1, 1), each=run), length(J.list))

Min.Alters <- bp_list$Min.Alters
Max.Alters <- bp_list$Max.Alters
Mean.Alters <- bp_list$Mean.Alters
Median.Alters <- bp_list$Median.Alters

outputs <- list(
  AOpts_one_trt_and_trt_fr_mis <- AOpts.df,
  YOpts_one_trt_and_trt_fr_mis <- YOpts.df,
  Regret_one_trt_and_trt_fr_mis <- Regret.df,
  Psi_Dw_one_trt_and_trt_fr_mis <- Psi.Dw,
  Psi_M1_one_trt_and_trt_fr_mis <- Psi.M1,
  Psi_M2_one_trt_and_trt_fr_mis <- Psi.M2,
  Psi_M3_one_trt_and_trt_fr_mis <- Psi.M3,
  Xsi_Dw_one_trt_and_trt_fr_mis <- Xsi.Dw,
  Xsi_M1_one_trt_and_trt_fr_mis <- Xsi.M1,
  Xsi_M2_one_trt_and_trt_fr_mis <- Xsi.M2,
  Xsi_M3_one_trt_and_trt_fr_mis <- Xsi.M3,
  Min_alters_one_trt_and_trt_fr_mis <- Min.Alters,
  Max_alters_one_trt_and_trt_fr_mis <- Max.Alters,
  Mean_alters_one_trt_and_trt_fr_mis <- Mean.Alters,
  Median_alters_one_trt_and_trt_fr_mis <- Median.Alters
)

saveRDS(outputs, "results/one_stage_trt_and_trt_fr_mis.rds")

