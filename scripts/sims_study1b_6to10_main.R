source(here::here("scripts", "sims_study1b_6to10.R")) 

run <- 1000
n <- 25
sam.size <- 25
a <- 2
b <- 2
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

bp_list <- begin.Sims(run, n=25, sam.size = 25, a=2, b=2, c=10, tau3=1, tau2=1, sigma2=1, alpha.0=-1, alpha.1=1, alpha.2=0.5, psi.0=1, psi.1=1, psi.2=-1, psi.3=-1.5)

J.list <- seq(a, b, 1)

phys.grps <- c()

for(j in J.list){
  phys.val <- as.character(j)
  #Times eleven because we have what is observed plus ten methods to compare
  phys.grps <- c(phys.grps, rep(phys.val, run*11))
}

AOpts.df <- bp_list$AOpts
AOpts.df$phys.grps <- phys.grps

YOpts.df <- bp_list$YOpts
YOpts.df$phys.grps <- phys.grps

YAgree.df <- bp_list$YAgree
YAgree.df$phys.grps <- phys.grps

phys.dw <- c()

for(j in J.list){
  phys.val <- as.character(j)
  #Times 4 because we have four blip parameters
  phys.dw <- c(phys.dw, rep(phys.val, run*4))
}

phys.dw.xsi <- c()

for(j in J.list){
    phys.val <- as.character(j)
    #Times 2 because we have two spillover parameters
    phys.dw.xsi <- c(phys.dw.xsi, rep(phys.val, run*2))
}

Psi.Dw <- bp_list$Psi.Dw
Psi.Dw$phys.grps <- phys.dw
Psi.Dw$TruePsi <- rep(rep(c(1, 1, -1, -1.5), each=run), length(J.list))

Psi.M1 <- bp_list$Psi.M1
Psi.M1$phys.grps <- phys.dw
#"Each" is based off simulation runs. Needs to be entered manually
Psi.M1$TruePsi <- rep(rep(c(1, 1, -1, -1.5), each=run), length(J.list))

Psi.M2 <- bp_list$Psi.M2
Psi.M2$phys.grps <- phys.dw
#"Each" is based off simulation runs. Needs to be entered manually
Psi.M2$TruePsi <- rep(rep(c(1, 1, -1, -1.5), each=run), length(J.list))

Psi.M3 <- bp_list$Psi.M3
Psi.M3$phys.grps <- phys.dw
Psi.M3$TruePsi <- rep(rep(c(1, 1, -1, -1.5), each=run), length(J.list))

Psi.M4 <- bp_list$Psi.M4
Psi.M4$phys.grps <- phys.dw
Psi.M4$TruePsi <- rep(rep(c(1, 1, -1, -1.5), each=run), length(J.list))

Psi.M5 <- bp_list$Psi.M5
Psi.M5$phys.grps <- phys.dw
Psi.M5$TruePsi <- rep(rep(c(1, 1, -1, -1.5), each=run), length(J.list))

Psi.M6 <- bp_list$Psi.M6
Psi.M6$phys.grps <- phys.dw
Psi.M6$TruePsi <- rep(rep(c(1, 1, -1, -1.5), each=run), length(J.list))

Psi.M7 <- bp_list$Psi.M7
Psi.M7$phys.grps <- phys.dw
Psi.M7$TruePsi <- rep(rep(c(1, 1, -1, -1.5), each=run), length(J.list))

Psi.M8 <- bp_list$Psi.M8
Psi.M8$phys.grps <- phys.dw
Psi.M8$TruePsi <- rep(rep(c(1, 1, -1, -1.5), each=run), length(J.list))

Psi.M9 <- bp_list$Psi.M9
Psi.M9$phys.grps <- phys.dw
Psi.M9$TruePsi <- rep(rep(c(1, 1, -1, -1.5), each=run), length(J.list))

Xsi.Dw <- bp_list$Xsi.Dw
Xsi.Dw$phys.grps <- phys.dw.xsi
#"Each" is based off simulation runs. Needs to be entered manually
Xsi.Dw$TrueXsi <- rep(rep(c(1, 1), each=run), length(J.list))

Xsi.M1 <- bp_list$Xsi.M1
Xsi.M1$phys.grps <- phys.dw.xsi
#"Each" is based off simulation runs. Needs to be entered manually
Xsi.M1$TrueXsi <- rep(rep(c(1, 1), each=run), length(J.list))

Xsi.M2 <- bp_list$Xsi.M2
Xsi.M2$phys.grps <- phys.dw.xsi
#"Each" is based off simulation runs. Needs to be entered manually
Xsi.M2$TrueXsi <- rep(rep(c(1, 1), each=run), length(J.list))

Xsi.M3 <- bp_list$Xsi.M3
Xsi.M3$phys.grps <- phys.dw.xsi
Xsi.M3$TrueXsi <- rep(rep(c(1, 1), each=run), length(J.list))

Xsi.M4 <- bp_list$Xsi.M4
Xsi.M4$phys.grps <- phys.dw.xsi
Xsi.M4$TrueXsi <- rep(rep(c(1, 1), each=run), length(J.list))

Xsi.M5 <- bp_list$Xsi.M5
Xsi.M5$phys.grps <- phys.dw.xsi
Xsi.M5$TrueXsi <- rep(rep(c(1, 1), each=run), length(J.list))

Xsi.M6 <- bp_list$Xsi.M6
Xsi.M6$phys.grps <- phys.dw.xsi
Xsi.M6$TrueXsi <- rep(rep(c(1, 1), each=run), length(J.list))

Xsi.M7 <- bp_list$Xsi.M7
Xsi.M7$phys.grps <- phys.dw.xsi
Xsi.M7$TrueXsi <- rep(rep(c(1, 1), each=run), length(J.list))

Xsi.M8 <- bp_list$Xsi.M8
Xsi.M8$phys.grps <- phys.dw.xsi
Xsi.M8$TrueXsi <- rep(rep(c(1, 1), each=run), length(J.list))

Xsi.M9 <- bp_list$Xsi.M9
Xsi.M9$phys.grps <- phys.dw.xsi
Xsi.M9$TrueXsi <- rep(rep(c(1, 1), each=run), length(J.list))

Mean.W.Alters <- bp_list$Mean.W.Alters
Mean.B.Alters <- bp_list$Mean.B.Alters

outputs <- list(
  AOpts_study1b_6to10 <- AOpts.df,
  YOpts_study1b_6to10 <- YOpts.df,
  Psi_Dw_study1b_6to10 <- Psi.Dw,
  Psi_M1_study1b_6to10 <- Psi.M1,
  Psi_M2_study1b_6to10 <- Psi.M2,
  Psi_M3_study1b_6to10 <- Psi.M3,
  Psi_M4_study1b_6to10 <- Psi.M4,
  Psi_M5_study1b_6to10 <- Psi.M5,
  Psi_M6_study1b_6to10 <- Psi.M6,
  Psi_M7_study1b_6to10 <- Psi.M7,
  Psi_M8_study1b_6to10 <- Psi.M8,
  Psi_M9_study1b_6to10 <- Psi.M9,

  Xsi_Dw_study1b_6to10 <- Xsi.Dw,
  Xsi_M1_study1b_6to10 <- Xsi.M1,
  Xsi_M2_study1b_6to10 <- Xsi.M2,
  Xsi_M3_study1b_6to10 <- Xsi.M3,
  Xsi_M4_study1b_6to10 <- Xsi.M4,
  Xsi_M5_study1b_6to10 <- Xsi.M5,
  Xsi_M6_study1b_6to10 <- Xsi.M6,
  Xsi_M7_study1b_6to10 <- Xsi.M7,
  Xsi_M8_study1b_6to10 <- Xsi.M8,
  Xsi_M9_study1b_6to10 <- Xsi.M9,
  
  Mean_W_alters_study1b_6to10 <- Mean.W.Alters,
  Mean_B_alters_study1b_6to10 <- Mean.B.Alters
)

saveRDS(outputs, "results/study1b_6to10.rds")

bp.trt <- ggplot(AOpts.df, aes(x=phys.grps, y=TrtRate, fill=Method)) + geom_boxplot()
bp.out <- ggplot(YOpts.df, aes(x=phys.grps, y=Outcome, fill=Method)) + geom_boxplot()
bp.M1 <- ggplot(Psi.M1, aes(x=PsiParameter, y=Estimate, fill=phys.grps)) + geom_boxplot()
bp.M2 <- ggplot(Psi.M2, aes(x=PsiParameter, y=Estimate, fill=phys.grps)) + geom_boxplot()
bp.M3 <- ggplot(Psi.M3, aes(x=PsiParameter, y=Estimate, fill=phys.grps)) + geom_boxplot()

bp.trt + labs(y="Optimal treatment rate", x="No. of physician groups") + scale_fill_brewer(palette = "Greys") + theme_classic()
bp.out + labs(y="Optimal outcome", x="No. of physician groups") + scale_fill_brewer(palette = "Greys") + theme_classic()

bp.M1 + labs(y="Estimate under ME2", x="Blip parameter") + scale_fill_brewer(palette = "Greys") + theme_classic()
bp.M2 + labs(y="Estimate under ME3", x="Blip parameter") + scale_fill_brewer(palette = "Greys") + theme_classic()
bp.M3 + labs(y="Estimate under ME4", x="Blip parameter") + scale_fill_brewer(palette = "Greys") + theme_classic()

bp.M1 + labs(y="Estimate under ME2", x="Spillover parameter") + scale_fill_brewer(palette = "Greys") + theme_classic()
bp.M2 + labs(y="Estimate under ME3", x="Spillover parameter") + scale_fill_brewer(palette = "Greys") + theme_classic()
bp.M3 + labs(y="Estimate under ME4", x="Spillover parameter") + scale_fill_brewer(palette = "Greys") + theme_classic()

plotdata1 <- na.omit(AOpts.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(TrtRate),
         mean = mean(TrtRate),
         sd = sd(TrtRate),
         empse = sqrt((sum((TrtRate-mean(TrtRate))^2))/(n-1)),
         .groups="keep")

plotdata1

plotdata2 <- na.omit(YOpts.df) %>%
  group_by(Method, phys.grps) %>%
  summarize(n = n(),
         median = median(Outcome),
         mean = mean(Outcome),
         sd = sd(Outcome),
         empse = sqrt((sum((Outcome-mean(Outcome))^2))/(n-1)),
         .groups="keep")

plotdata2

plotdata3 <- na.omit(Psi.M2) %>%
  group_by(PsiParameter, phys.grps) %>%
  summarize(n = n(),
         median = median(Estimate),
         mean = mean(Estimate),
         empse = sqrt((sum((Estimate-mean(Estimate))^2))/(n-1)),
         mse = (sum((Estimate-TruePsi)^2))/n,
         .groups="keep")

plotdata3

plotdata4 <- na.omit(Psi.M3) %>%
  group_by(PsiParameter, phys.grps) %>%
  summarize(n = n(),
         median = median(Estimate),
         mean = mean(Estimate),
         empse = sqrt((sum((Estimate-mean(Estimate))^2))/(n-1)),
         mse = (sum((Estimate-TruePsi)^2))/n,
         .groups="keep")

plotdata4

plotdata5 <- na.omit(Xsi.M1) %>%
  group_by(XsiParameter, phys.grps) %>%
  summarize(n = n(),
         median = median(Estimate),
         mean = mean(Estimate),
         empse = sqrt((sum((Estimate-mean(Estimate))^2))/(n-1)),
         mse = (sum((Estimate-TrueXsi)^2))/n,
         .groups="keep")

plotdata5

plotdata6 <- na.omit(Xsi.M2) %>%
  group_by(XsiParameter, phys.grps) %>%
  summarize(n = n(),
         median = median(Estimate),
         mean = mean(Estimate),
         empse = sqrt((sum((Estimate-mean(Estimate))^2))/(n-1)),
         mse = (sum((Estimate-TrueXsi)^2))/n,
         .groups="keep")

plotdata6

plotdata7 <- na.omit(Xsi.M3) %>%
  group_by(XsiParameter, phys.grps) %>%
  summarize(n = n(),
         median = median(Estimate),
         mean = mean(Estimate),
         empse = sqrt((sum((Estimate-mean(Estimate))^2))/(n-1)),
         mse = (sum((Estimate-TrueXsi)^2))/n,
         .groups="keep")

plotdata7

