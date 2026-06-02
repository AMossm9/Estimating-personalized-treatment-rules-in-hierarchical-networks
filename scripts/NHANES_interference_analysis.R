#install.packages("nhanesA")
library(nhanesA)
library(dplyr)
library(lme4)
library(ggplot2)
library(DTRreg)
#_J represents the 2017-2018 cycle from NHANES
#DEMO_J represents demographics
#ALQ_J represents alcohol consumption questionnaire
#Outcome of interest in ALQ_J: no. of days have 4 or 5 drinks/past 12 mos
#BPQ_J represents awareness, treatment, and control of high blood pressure
#and cholesterol
#HOQ_J represents household characteristics

data_gen <- function(seed, df_cleaned){
  
  #Generate hierarchical data
  I <- 2
  J <- 4
  L3 <- rep(c(1, 2), c(172, 263))
  L2 <- rep(c(1, 2, 3, 4, 5, 6, 7, 8), c(43, 43, 43, 43, 65, 66, 66, 66))
  i_counts <- as.vector(table(L3))
  
  j_counts <- as.vector(table(L2))
  
  b_0i <- rnorm(I, 0, 1) 
  
  b_0j <- rnorm(I*J, 0, 1) 
  
  beta_0i <- rep(b_0i, i_counts)
  
  beta_0j <- rep(b_0j, j_counts)
  
  grouping_cols <- as.data.frame(cbind(L3, L2, beta_0i, beta_0j))
  
  df_cleaned_hous <- df_cleaned %>% arrange(HOUSING)
  
  df_nested <- as.data.frame(cbind(grouping_cols, df_cleaned_hous))
  
  df_nested$HOUSING <- ifelse(df_nested$HOUSING=="Rented", 1, 0)
  
  df_nested$GENDER <- ifelse(df_nested$GENDER=="Male", 1, 0)
  
  df_nested$EDUC <- ifelse(df_nested$EDUC=="College graduate or above" | df_nested$EDUC=="Some college or AA degree", 1, 0)
  
  df_nested <- df_nested[!is.na(df_nested$EDUC),]
  
  logit.p.1 <- 4 - 0.20*df_nested$AGE + 2*df_nested$HOUSING + df_nested$EDUC + 0.05*df_nested$GENDER + df_nested$beta_0i + df_nested$beta_0j
  
  logit.p.2 <- 4 - 0.20*df_nested$AGE + 2*df_nested$HOUSING + df_nested$EDUC + 0.05*df_nested$GENDER + df_nested$beta_0i + df_nested$beta_0j
  
  prob.trt.1 <- exp(logit.p.1)/(exp(logit.p.1)+1)
  A_1 <- rbinom(nrow(df_nested), 1, prob.trt.1)
  
  prob.trt.2 <- exp(logit.p.2)/(exp(logit.p.2)+1)
  A_2 <- rbinom(nrow(df_nested), 1, prob.trt.2)
  
  trt.indicator.1 <- ifelse(A_1 == 1, 1, -1)
  trt.indicator.2 <- ifelse(A_2 == 1, 1, -1)
  
  df_nested$A1 <- A_1
  df_nested$A2 <- A_2
  
  df_nested$trt.indicator.1 <- ifelse(A_1 == 1, 1, -1)
  df_nested$trt.indicator.2 <- ifelse(A_2 == 1, 1, -1)
  
  
  #Fixed- and random-effects model, stage 1
  mod.1 <- glmer(A1~AGE + GENDER + HOUSING + EDUC + (1 | L3/L2), family = 'binomial', data = df_nested)
  
  #Fixed- and random-effects model, stage 2
  mod.2 <- glmer(A2~AGE + GENDER + HOUSING + EDUC + (1 | L3/L2), family = 'binomial', data = df_nested)
  
  df_nested$logit.p.1 <- predict(mod.1, df_nested, type='response')
  
  df_nested$logit.p.2 <- predict(mod.2, df_nested, type='response')
  
  return(df_nested)
}

egos.alters <- function(seed, df_nested, I, J, n){
  
  trt.data <- df_nested
  
  set.seed(seed)
  
  sam.size <- 50
  
  mydata.egos <- list()
  mydata.alters <- list()
  
  #Subset data based on Level 3
  for(i in 1:I){
    sam.size <- ifelse(i == 1, 86, 131)
    df.L3 <- trt.data[trt.data$L3 == i,]
    
    #Total no. of rows = total no. of unique unit ID's
    #Randomly sample egos from available ID's, the number of which 
    #must be less than or equal to half of total units in this Level 3 group
    #This is part of ensuring that no ego is left without an alter
    IDs <- df.L3$ID
    egos <- sample(IDs, sam.size, replace=FALSE)
    
    #The remaining individuals in this Level 3 hierarchy become alters
    alters <- IDs[! IDs %in% egos]
    
    #Randomly assign alters to egos so that each ego has at least one alter
    set.seed(seed)
    X <- split(alters, c(seq(1:length(egos)), sample(length(egos), length(alters)-length(egos), repl=TRUE)))
    #print(X)
    id.egos <- c()
    id.alters <- c()
    egodatalist <- list()
    alterdatalist <- list()
    
    #Select the rows corresponding to egos in the sample, and repeat each of them
    #based on how many alters that ego has. Collect these as a list. 
    #Collect the rows corresponding to alters for each ego as a list.
    #Finally, stack the egos' rows, stack the alters' rows, and merge them using cbind().
    #This way, we have rows for each ego-alter assignment.
    for(c in 1:length(egos)){
      e <- egos[c]
      
      id <- rep(e, length(X[[c]]))
      
      id.egos <- c(id.egos, id)
      
      id.alters <- c(id.alters, X[[c]])
      
      df.egos.sam <- do.call(rbind, replicate(length(X[[c]]), df.L3[df.L3$ID==e,], simplify = FALSE))
      
      egodatalist[[c]] <- df.egos.sam
      
      alterdatalist[[c]] <- df.L3[df.L3$ID %in% X[[c]],]
    }
    
    df.egos.rep <- do.call(rbind, egodatalist)
    df.alters <- do.call(rbind, alterdatalist)
    df.egos.rep$ego <- rep("ego", nrow(df.egos.rep))
    #names(df.egos.rep)[2] <- "ego.L2"
    df.alters$alter <- rep("alter", nrow(df.alters))
    #names(df.alters)[2] <- "alter.L2"
    
    mydata.egos[[i]] <- df.egos.rep
    #print(mydata.egos[[i]])
    mydata.alters[[i]] <- df.alters
    #print(mydata.alters[[i]])
  }
  
  data.egos <- do.call(rbind, mydata.egos)
  names(data.egos)[2] <- "ego.L2"
  data.alters <- do.call(rbind, mydata.alters)
  names(data.alters)[2] <- "alter.L2"
  names(data.alters)[-c(1,2)] <- paste(names(data.alters)[-c(1,2)], ".alter", sep="")
  df.network <- cbind(data.egos, data.alters)
  
  #Count the number of within- and between-group neighbours for each ego
  #based on the L2 group that each alter belongs to.
  #print(df.network)
  wthn.btwn.group <- c()
  
  for(i in 1:nrow(df.network)){
    if(df.network$ego.L2[i] == df.network$alter.L2[i]){
      wthn.btwn.group <- c(wthn.btwn.group, "W")
      
    }
    else{
      wthn.btwn.group <- c(wthn.btwn.group, "B")
      
    }
  }
  
  df.network$alter.group <- wthn.btwn.group
  mydata <- df.network
  
  #Should return a dataframe that has 25 columns
  return(mydata)
}


weight_calc <- function(seed, mydata, I, J, n){
  
  trt.data <- mydata
  
  alters.joint.1 <- c()
  total.w.1 <- c()
  total.b.1 <- c()
  
  alters.joint.2 <- c()
  total.w.2 <- c()
  total.b.2 <- c()
  
  #Iterate through egos, first through Level 3 and then through Level 2 in Level 3
  for(i in 1:I){
    
    df.L3 <- trt.data[trt.data$L3 == i,]
    ego.w.1 <- c()
    ego.b.1 <- c()
    
    ego.w.2 <- c()
    ego.b.2 <- c()
    
    #We need to search through each unique ego ID in our dataframe (since
    #there are multiple rows for an ego, each linked with an alters'information).
    unique.egos <- unique(df.L3$ID)
    for(u in unique.egos){
      
      #Extract dataframe for that ego
      df.ego <- df.L3[df.L3$ID==u,]
      
      #Find its within-group neighbours
      df.ego.w <- df.ego[df.ego$alter.group == "W",]
      
      sum.trt.w.1 <- sum(df.ego.w$A1.alter)
      
      sum.trt.w.2 <- sum(df.ego.w$A2.alter)
      
      #Find its between-group neighbours
      df.ego.b <- df.ego[df.ego$alter.group == "B",]
      
      sum.trt.b.1 <- sum(df.ego.b$A1.alter)
      
      sum.trt.b.2 <- sum(df.ego.b$A2.alter)
      
      ego.w.1 <- c(ego.w.1, rep(sum.trt.w.1, nrow(df.ego)))
      
      ego.b.1 <- c(ego.b.1, rep(sum.trt.b.1, nrow(df.ego)))
      
      ego.w.2 <- c(ego.w.2, rep(sum.trt.w.2, nrow(df.ego)))
      
      ego.b.2 <- c(ego.b.2, rep(sum.trt.b.2, nrow(df.ego)))
      
      
    }
    
    total.w.1 <- c(total.w.1, ego.w.1)
    total.b.1 <- c(total.b.1, ego.b.1)
    total.w.2 <- c(total.w.2, ego.w.2)
    total.b.2 <- c(total.b.2, ego.b.2)
    
  }
  
  trt.data$total.w.1 <- total.w.1
  trt.data$total.b.1 <- total.b.1
  trt.data$total.w.2 <- total.w.2
  trt.data$total.b.2 <- total.b.2
  #This is where alters.joint.1 and alters.joint.2 turn into the total number of treated neighbours
  trt.data$total.t.1 <- total.w.1 + total.b.1
  trt.data$total.t.2 <- total.w.2 + total.b.2
  
  #w.L3: network absolute value weights that account for both within- and between-group neighbours
  #w.within: network absolute value weights that account for only within-group neighbours
  #w.standard: calculates the ego's absolute value weight using fixed-effects only model,
  #does not account for neighbours
  #w.within: network absolute value weights that account for only within-group neighbours
  
  w.L3.1 <- c()
  w.within.1 <- c()
  w.standard.1 <- c()
  w.network.1 <- c()
  
  w.L3.2 <- c()
  w.within.2 <- c()
  w.standard.2 <- c()
  w.network.2 <- c()
  
  #Subset data based on Level 3
  for(i in 1:I){
    df.L3 <- trt.data[trt.data$L3 == i,]
    #Use ".L2" as we start storing weights
    #because we are now working with 
    #egos in Level 2 groups within a given
    #Level 3 group
    w.L2.1 <- c()
    w.within.L2.1 <- c()
    w.standard.L2.1 <- c()
    w.network.L2.1 <- c()
    
    w.L2.2 <- c()
    w.within.L2.2 <- c()
    w.standard.L2.2 <- c()
    w.network.L2.2 <- c()
    
    unique.egos <- unique(df.L3$ID)
    
    #Iterate through each ego within Level 3 (possibly
    #from different Level 2 groups)
    for(u in unique.egos){
      
      #Get dataframe corresponding to ego
      df.ego <- df.L3[df.L3$ID==u,]
      
      #Get the ego's probability of treatment from stage 1 mixed model (it's called logit but is actually the expit)
      #and indicator for whether or not they received treatment "1"
      ego.logit.1 <- df.ego$logit.p.1[1]
      ego.indicator.1 <- df.ego$trt.indicator.1[1]
      ego.A.1 <- ifelse(ego.indicator.1==1, 1, 0)
      
      #Get vector of alters' (within- and between-group) probability of treatment from stage 1 mixed model
      alter.logit.1 <- df.ego$logit.p.1.alter
      alter.indicator.1 <- df.ego$trt.indicator.1.alter
      alter.A.1 <- ifelse(alter.indicator.1==1, 1, 0)
      
      #Get the ego's probability of treatment from stage 2 mixed model (it's called logit but is actually the expit)
      #and indicator for whether or not they received treatment "1"
      ego.logit.2 <- df.ego$logit.p.2[1]
      ego.indicator.2 <- df.ego$trt.indicator.2[1]
      ego.A.2 <- ifelse(ego.indicator.2==1, 1, 0)
      
      #Get vector of alters' (within- and between-group) probability of treatment from stage 2 mixed model
      alter.logit.2 <- df.ego$logit.p.2.alter
      alter.indicator.2 <- df.ego$trt.indicator.2.alter
      alter.A.2 <- ifelse(alter.indicator.2==1, 1, 0)
      
      #Calculate the stage 1 and stage 2 network absolute value Weights for group propensity score
      w.L2.1 <- c(w.L2.1, rep(abs(ego.A.1-ego.logit.1)*prod(abs(alter.A.1-alter.logit.1)), nrow(df.ego)))
      w.L2.2 <- c(w.L2.2, rep(abs(ego.A.2-ego.logit.2)*prod(abs(alter.A.2-alter.logit.2)), nrow(df.ego)))
      
      #IMPORTANT:
      #Refer to these links and the formula (3) in Perez-Heydrich et al. (2014) Paper to understand why we integrated out random effects:
      #using re.form=NA
      #https://stackoverflow.com/questions/36153245/non-finite-function-value-error-in-integrate-function-in-r
      #https://stats.stackexchange.com/questions/308927/re-form-specifying-one-random-effect-but-not-the-other
      
      #Weights that will account for within-group neighbours only
      #Check number of within-group neighbours
      #If > 0, take product of ego and within-group alters' propensities from model with mixed effects
      #Else, it's just the absolute value weight of the ego's propensity (but still from the mixed effects model)
      df.ego.w <- df.ego[df.ego$alter.group == "W",]
      if(nrow(df.ego.w)==0){
        ego.w.logit.1 <- df.ego$logit.p.1[1]
        ego.w.indicator.1 <- df.ego$trt.indicator.1[1]
        ego.w.A.1 <- ifelse(ego.w.indicator.1==1, 1, 0)
        
        ego.w.logit.2 <- df.ego$logit.p.2[1]
        ego.w.indicator.2 <- df.ego$trt.indicator.2[1]
        ego.w.A.2 <- ifelse(ego.w.indicator.2==1, 1, 0)
        
        w.within.L2.1 <- c(w.within.L2.1, rep(abs(ego.w.A.1-ego.w.logit.1), nrow(df.ego)))
        w.network.L2.1 <- c(w.network.L2.1, rep(abs(ego.w.A.1-ego.w.logit.1), nrow(df.ego)))
        
        w.within.L2.2 <- c(w.within.L2.2, rep(abs(ego.w.A.2-ego.w.logit.2), nrow(df.ego)))
        w.network.L2.2 <- c(w.network.L2.2, rep(abs(ego.w.A.2-ego.w.logit.2), nrow(df.ego)))
        
      }
      
      else{
        ego.w.logit.1 <- df.ego$logit.p.1[1]
        ego.w.indicator.1 <- df.ego$trt.indicator.1[1]
        ego.w.A.1 <- ifelse(ego.w.indicator.1==1, 1, 0)
        
        alter.w.logit.1 <- df.ego.w$logit.p.1.alter
        alter.w.indicator.1 <- df.ego.w$trt.indicator.1.alter
        alter.w.A.1 <- ifelse(alter.w.indicator.1==1, 1, 0)
        
        ego.w.logit.2 <- df.ego$logit.p.2[1]
        ego.w.indicator.2 <- df.ego$trt.indicator.2[1]
        ego.w.A.2 <- ifelse(ego.w.indicator.2==1, 1, 0)
        
        alter.w.logit.2 <- df.ego.w$logit.p.2.alter
        alter.w.indicator.2 <- df.ego.w$trt.indicator.2.alter
        alter.w.A.2 <- ifelse(alter.w.indicator.2==1, 1, 0)
        
        w.within.L2.1 <- c(w.within.L2.1, rep(abs(ego.w.A.1-ego.w.logit.1)*prod(abs(alter.w.A.1-alter.w.logit.1)), nrow(df.ego)))
        w.network.L2.1 <- c(w.network.L2.1, rep(abs(ego.w.A.1-ego.w.logit.1)*prod(abs(alter.w.A.1-alter.w.logit.1)), nrow(df.ego)))
        
        w.within.L2.2 <- c(w.within.L2.2, rep(abs(ego.w.A.2-ego.w.logit.2)*prod(abs(alter.w.A.2-alter.w.logit.2)), nrow(df.ego)))
        w.network.L2.2 <- c(w.network.L2.2, rep(abs(ego.w.A.2-ego.w.logit.2)*prod(abs(alter.w.A.2-alter.w.logit.2)), nrow(df.ego)))
        
      }
      
      #Weights for standard individual absolute value weight
      #from mixed-effects model
      ego.s.logit.1 <- df.ego$logit.p.1[1]
      ego.s.indicator.1 <- df.ego$trt.indicator.1[1]
      ego.s.A.1 <- ifelse(ego.s.indicator.1==1, 1, 0)
      
      ipw.s.1 <- abs(ego.s.A.1 - ego.s.logit.1)
      w.standard.L2.1 <- c(w.standard.L2.1, rep(ipw.s.1, nrow(df.ego)))
      
      ego.s.logit.2 <- df.ego$logit.p.2[1]
      ego.s.indicator.2 <- df.ego$trt.indicator.2[1]
      ego.s.A.2 <- ifelse(ego.s.indicator.2==1, 1, 0)
      ipw.s.2 <- abs(ego.s.A.2 - ego.s.logit.2)
      w.standard.L2.2 <- c(w.standard.L2.2, rep(ipw.s.2, nrow(df.ego)))
    }
    
    w.L3.1 <- c(w.L3.1, w.L2.1)
    w.within.1 <- c(w.within.1, w.within.L2.1)
    w.standard.1 <- c(w.standard.1, w.standard.L2.1)
    w.network.1 <- c(w.network.1, w.network.L2.1)
    
    w.L3.2 <- c(w.L3.2, w.L2.2)
    w.within.2 <- c(w.within.2, w.within.L2.2)
    w.standard.2 <- c(w.standard.2, w.standard.L2.2)
    w.network.2 <- c(w.network.2, w.network.L2.2)
  }
  
  trt.data$w.L3.1 <- w.L3.1
  trt.data$w.within.1 <- w.within.1
  trt.data$w.standard.1 <- w.standard.1
  trt.data$w.network.1 <- w.network.1
  
  trt.data$w.L3.2 <- w.L3.2
  trt.data$w.within.2 <- w.within.2
  trt.data$w.standard.2 <- w.standard.2
  trt.data$w.network.2 <- w.network.2
  
  return(trt.data)
  
}

main_analysis <- function(seed, df_cleaned, I, J, n){
  
  df_nested <- data_gen(seed, df_cleaned)
  
  mydata <- egos.alters(seed, df_nested, I, J, n)
  
  trt.data <- weight_calc(seed, mydata, I, J, n)
  
  ego.df <- trt.data
  
  num.reps <- as.numeric(table(ego.df$ID))
  
  ego.df <- ego.df[order(ego.df$ID),]
  
  ego.df$num.alters <- rep(num.reps, num.reps) 
  
  ego.df <- ego.df[!duplicated(ego.df$ID), ]
  
  ego.df <- subset(ego.df, select=c(ID, L3:num.alters))
  
  ego.df <- ego.df[order(ego.df$num.alters),]
  
  return(list(num_original = nrow(mydata), trt.data=ego.df))
}

run_dwols_intf <- function(df_dwols){
  #Implement dWOLS using DTRreg package
  
  alc_freq <- df_dwols$alc_freq
  
  AGE <- df_dwols$AGE
  
  HOUSING <- df_dwols$HOUSING
  
  A1 <- df_dwols$A1
  
  A2 <- df_dwols$A2
  
  treat_wgt_man <- list(stage_1 = df_dwols$w.standard.1, stage_2 = df_dwols$w.standard.2)
  
  blip.mod <- list(~ total.t.1 + HOUSING, ~ total.t.2+HOUSING)
  
  treat.mod <- list(A1 ~ 1, A2 ~ 1)
  
  tf.mod <- list(~ total.t.1 + AGE + HOUSING, ~ total.t.2 + AGE + HOUSING)
  
  mod1 <- DTRreg(alc_freq, blip.mod, treat.mod, tf.mod,
                 data = df_dwols, method = "dwols", weight = "manual", treat.wgt.man = treat_wgt_man, var.estim="bootstrap")
  
  return(mod1)
}


run_dwols_intf_wthn <- function(df_dwols){
  #Implement dWOLS using DTRreg package
  
  alc_freq <- df_dwols$alc_freq
  
  AGE <- df_dwols$AGE
  
  HOUSING <- df_dwols$HOUSING
  
  A1 <- df_dwols$A1
  
  A2 <- df_dwols$A2
  
  treat_wgt_man <- list(stage_1 = df_dwols$w.standard.1, stage_2 = df_dwols$w.standard.2)
  
  blip.mod <- list(~ total.w.1 + HOUSING, ~ total.w.2+HOUSING)
  
  treat.mod <- list(A1 ~ 1, A2 ~ 1)
  
  tf.mod <- list(~ total.w.1 + AGE + HOUSING, ~ total.w.2 + AGE + HOUSING)
  
  mod1 <- DTRreg(alc_freq, blip.mod, treat.mod, tf.mod,
                 data = df_dwols, method = "dwols", weight = "manual", treat.wgt.man = treat_wgt_man, var.estim="bootstrap")
  
  return(mod1)
}


num_sims <- 100

num_patients_original_train <- matrix(NA, num_sims, 1)

num_patients_original_test <- matrix(NA, num_sims, 1)

num_patients_train <- matrix(NA, num_sims, 1)

num_patients_test <- matrix(NA, num_sims, 1)


AOpts_test <- matrix(NA, num_sims, 2)

blips <- matrix(NA, num_sims, 6)

spillovers <- matrix(NA, num_sims, 2)

regrets <- matrix(NA, num_sims, 3)

predicted_outcomes <- matrix(NA, num_sims, 2)


AOpts_age_stage1 <- matrix(NA, num_sims, 4)

AOpts_age_stage2 <- matrix(NA, num_sims, 4)

ANonOpts_age_stage1 <- matrix(NA, num_sims, 4)

ANonOpts_age_stage2 <- matrix(NA, num_sims, 4)


AOpts_test_wthn <- matrix(NA, num_sims, 2)

blips_wthn <- matrix(NA, num_sims, 6)

spillovers_wthn <- matrix(NA, num_sims, 2)

regrets_wthn <- matrix(NA, num_sims, 3)

predicted_outcomes_wthn <- matrix(NA, num_sims, 2)


AOpts_age_stage1_wthn <- matrix(NA, num_sims, 4)

AOpts_age_stage2_wthn <- matrix(NA, num_sims, 4)

ANonOpts_age_stage1_wthn <- matrix(NA, num_sims, 4)

ANonOpts_age_stage2_wthn <- matrix(NA, num_sims, 4)


#Running the simulations: first, clean the data from 2017-2018 cycle of NHANES
demo <- nhanes('DEMO_J')
alc <- nhanes('ALQ_J') #142 for no. of days have 4 or 5 drinks/past 12 mos
bp <- nhanes('BPQ_J') #BPQ050A for taking prescribed medication for HBP
hous <- nhanes('HOQ_J') #HOQ065 - Home owned, bought, rented, other

demo_translate <- nhanesTranslate('DEMO_J',
                                  c('SEQN', # Respondent sequence number
                                    'RIAGENDR', # Gender
                                    'RIDAGEYR', # Age in years at screening
                                    'RIDRETH3', # Race
                                    'DMDEDUC2'), # Education level - Adults 20+
                                  data = demo)

new_demo <- select(demo_translate,
                   c(SEQN,
                     RIAGENDR,
                     RIDAGEYR,
                     RIDRETH3,
                     DMDEDUC2))
head(new_demo,5)

hous_translate <- nhanesTranslate('HOQ_J',
                                  c('SEQN', 
                                    'HOD050', 
                                    'HOQ065'),
                                  data = hous)

new_hous <- select(hous_translate, c(SEQN, HOQ065))

bp_translate <- nhanesTranslate('BPQ_J',
                                c('SEQN', 
                                  'BPQ050A'),
                                data = bp)

new_bp <- select(bp_translate, c(SEQN, BPQ050A))


alc_translate <- nhanesTranslate('ALQ_J',
                                 c('SEQN', 
                                   'ALQ142'),
                                 data = alc)

new_alc <- select(alc_translate, c(SEQN, ALQ142))

demo_hous <- left_join(new_demo, new_hous, by = "SEQN")

demo_hous_bp <- left_join(demo_hous, new_bp, by = "SEQN")

demo_hous_bp_alc <- left_join(demo_hous_bp, new_alc, by = "SEQN")

full_df <- demo_hous_bp_alc

names(full_df) <- c("ID", "GENDER", "AGE", "ETHN", "EDUC", "HOUSING", "BP_MED", "ALC")

summary(full_df$AGE)

df_filt_age <- filter(full_df, AGE <= 28 & AGE >= 21)

df_filt_hous <- filter(df_filt_age, HOUSING == "Rented" | HOUSING == "Owned or being bought")

df_filt_alc <- filter(df_filt_hous, ALC != "Every day" & ALC != "Nearly every day" & ALC != "3 to 4 times a week" & ALC != "Don't know" & ALC != "Refused")

df_filt_alc_new <- df_filt_alc %>%
  mutate(
    alc_freq = case_when(
      ALC == 'Never in the last year' ~ 8,
      ALC == '1 to 2 times in the last year' ~ 7,
      ALC == '3 to 6 times in the last year' ~ 6,
      ALC == '7 to 11 times in the last year' ~ 5,
      ALC == 'Once a month' ~ 4,
      ALC == '2 to 3 times a month' ~ 3,
      ALC == 'Once a week' ~ 2,
      ALC == '2 times a week' ~ 1)
  )

df_cleaned <- df_filt_alc_new[!is.na(df_filt_alc_new$alc_freq) & !is.na(df_filt_alc_new$ALC), ]

I <- 2

J <- 4

n <- 50

for(i in 1:num_sims){

  seed <- i
  
  print(seed)
  
  df_dwols <- main_analysis(seed, df_cleaned, I, J, n)$trt.data
  
  num_patients_train[i,1] <- nrow(df_dwols)
  
  num_patients_original_train[i,1] <- main_analysis(seed, df_cleaned, I, J, n)$num_original
  
  df_test <- main_analysis(seed+1, df_cleaned, I, J, n)$trt.data
  
  num_patients_test[i,1] <- nrow(df_test)
  
  num_patients_original_test[i,1] <- main_analysis(seed+1, df_cleaned, I, J, n)$num_original
  
  dwols_res <- run_dwols_intf(df_dwols)
  
  blip_stage1 <- c(dwols_res$psi[[1]][1], dwols_res$psi[[1]][2], dwols_res$psi[[1]][3])
  
  print(blip_stage1)
  
  blip_stage2 <- c(dwols_res$psi[[2]][1], dwols_res$psi[[2]][2], dwols_res$psi[[2]][3])
  
  print(blip_stage2)
  
  xi_stage1 <- c(dwols_res$beta[[1]][2])
  
  xi_stage2 <- c(dwols_res$beta[[2]][2])
  
  AOpt_stage1 <- as.numeric((blip_stage1[[1]] + blip_stage1[[2]]*df_test$total.t.1 + blip_stage1[[3]]*df_test$HOUSING) > 0)
  
  AOpt_stage2 <- as.numeric((blip_stage2[1] + blip_stage2[2]*df_test$total.t.2 + blip_stage2[3]*df_test$HOUSING) > 0)
  
  regret_stage1 <- mean((AOpt_stage1-df_test$A1)*(blip_stage1[1] + blip_stage1[2]*df_test$total.t.1 + blip_stage1[3]*df_test$HOUSING))
  
  regret_stage2 <- mean((AOpt_stage2-df_test$A2)*(blip_stage2[1] + blip_stage2[2]*df_test$total.t.2 + blip_stage2[3]*df_test$HOUSING))
  
  regret_total <- mean((AOpt_stage1-df_test$A1)*(blip_stage1[1] + blip_stage1[2]*df_test$total.t.1 + blip_stage1[3]*df_test$HOUSING) + (AOpt_stage2-df_test$A2)*(blip_stage2[1] + blip_stage2[2]*df_test$total.t.2 + blip_stage2[3]*df_test$HOUSING))
  
  #Store results in matrices
  
  AOpts_test[i,1] <- mean(AOpt_stage1==df_test$A1)
  
  AOpts_test[i,2] <- mean(AOpt_stage2==df_test$A2)
  
  
  blips[i, 1] <- blip_stage1[1]
  
  blips[i, 2] <- blip_stage1[2]
  
  blips[i, 3] <- blip_stage1[3]
  
  blips[i, 4] <- blip_stage2[1]
  
  blips[i, 5] <- blip_stage2[2]
  
  blips[i, 6] <- blip_stage2[3]
  
  spillovers[i, 1] <- xi_stage1[1]
  
  spillovers[i, 2] <- xi_stage2[1]
  
  
  regrets[i,1] <- regret_stage1
  
  regrets[i,2] <- regret_stage2
  
  regrets[i,3] <- regret_total
  
  #predicted_outcomes[i,1] <- mean(predict(dwols_res, df_test))
  
  #predicted_outcomes[i,2] <- median(predict(dwols_res, df_test))
  
  dwols_res_wthn <- run_dwols_intf_wthn(df_dwols)
  
  blip_stage1_wthn <- c(dwols_res_wthn$psi[[1]][1], dwols_res_wthn$psi[[1]][2], dwols_res_wthn$psi[[1]][3])
  
  print(blip_stage1_wthn)
  
  blip_stage2_wthn <- c(dwols_res_wthn$psi[[2]][1], dwols_res_wthn$psi[[2]][2], dwols_res_wthn$psi[[2]][3])
  
  print(blip_stage2_wthn)
  
  xi_stage1_wthn <- c(dwols_res_wthn$beta[[1]][2])
  
  xi_stage2_wthn <- c(dwols_res_wthn$beta[[2]][2])
  
  
  AOpt_stage1_wthn <- as.numeric((blip_stage1_wthn[1] + blip_stage1_wthn[2]*df_test$total.w.1 + blip_stage1_wthn[3]*df_test$HOUSING) > 0)
  
  AOpt_stage2_wthn <- as.numeric((blip_stage2_wthn[1] + blip_stage2_wthn[2]*df_test$total.w.2 + blip_stage2_wthn[3]*df_test$HOUSING) > 0)
  
  regret_stage1_wthn <- mean((AOpt_stage1_wthn-df_test$A1)*(blip_stage1_wthn[1] + blip_stage1_wthn[2]*df_test$total.w.1 + blip_stage1_wthn[3]*df_test$HOUSING))
  
  regret_stage2_wthn <- mean((AOpt_stage2_wthn-df_test$A2)*(blip_stage2_wthn[1] + blip_stage2_wthn[2]*df_test$total.w.2 + blip_stage2_wthn[3]*df_test$HOUSING))
  
  regret_total_wthn <- mean((AOpt_stage1_wthn-df_test$A1)*(blip_stage1_wthn[1] + blip_stage1_wthn[2]*df_test$total.w.1 + blip_stage1_wthn[3]*df_test$HOUSING) + (AOpt_stage2_wthn-df_test$A2)*(blip_stage2_wthn[1] + blip_stage2_wthn[2]*df_test$total.w.2 + blip_stage2_wthn[3]*df_test$HOUSING))
  
  #Store results in matrices
  
  AOpts_test_wthn[i,1] <- mean(AOpt_stage1_wthn==df_test$A1)
  
  AOpts_test_wthn[i,2] <- mean(AOpt_stage2_wthn==df_test$A2)
  
  
  blips_wthn[i, 1] <- blip_stage1_wthn[1]
  
  blips_wthn[i, 2] <- blip_stage1_wthn[2]
  
  blips_wthn[i, 3] <- blip_stage1_wthn[3]
  
  blips_wthn[i, 4] <- blip_stage2_wthn[1]
  
  blips_wthn[i, 5] <- blip_stage2_wthn[2]
  
  blips_wthn[i, 6] <- blip_stage2_wthn[3]
  
  spillovers_wthn[i, 1] <- xi_stage1_wthn[1]
  
  spillovers_wthn[i, 2] <- xi_stage2_wthn[1]
  
  
  regrets_wthn[i,1] <- regret_stage1_wthn
  
  regrets_wthn[i,2] <- regret_stage2_wthn
  
  regrets_wthn[i,3] <- regret_total_wthn
  
  #predicted_outcomes_wthn[i,1] <- mean(predict(dwols_res_wthn, df_test))
  
  #predicted_outcomes_wthn[i,2] <- median(predict(dwols_res_wthn, df_test))
}

outputs_btwn_and_wthn <- list(
  nhanes_intf_btwn_and_wthn_num_patients_train <- num_patients_train,
  nhanes_intf_btwn_and_wthn_num_patients_test <- num_patients_test,

  nhanes_intf_btwn_and_wthn_num_patients_original_train <- num_patients_original_train,
  nhanes_intf_btwn_and_wthn_num_patients_original_test <- num_patients_original_test,

  nhanes_intf_btwn_and_wthn_AOpts_test <- AOpts_test,
  nhanes_intf_btwn_and_wthn_blips <- blips,

  nhanes_intf_btwn_and_wthn_spillovers <- spillovers,

  nhanes_intf_btwn_and_wthn_regrets <- regrets
)

saveRDS(outputs_btwn_and_wthn, "results/nhanes_intf_btwn_and_wthn")

outputs_wthn <- list(

  nhanes_intf_wthn_AOpts_test <- AOpts_test,
  nhanes_intf_wthn_blips <- blips,

  nhanes_intf_wthn_spillovers <- spillovers,

  nhanes_intf_wthn_regrets <- regrets
)

saveRDS(outputs_wthn, "results/nhanes_intf_wthn")
#*************************************************************************************************************************************

