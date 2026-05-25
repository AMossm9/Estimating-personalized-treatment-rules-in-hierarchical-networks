#**********DWOLS FOR WITHIN- AND BETWEEN-GROUP INTERFERENCE: Two-STAGE SETTING**********

#This script simulates multiple hierarchies of network data, where egos are randomly selected from
#each hierarchy and randomly assigned neighbours within- and between-subgroups in the hierarchy.

#Install packages for mixed effects models, data visualization with ggplot,
#and incorporating latex expressions in figures

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

#**********DATA GENERATION**********

#Generate hierarchical, nested data
#Level 1 (the units or individuals) is nested in Level 2 (e.g. a physician group, household, etc.) which is 
#nested in Level 3 (e.g. a clinic, larger community, etc.)

#I = no. of Level 3 units (i.e. hierarchies)
#J = no. of Level 2 units PER Level 3 unit (meaning there are J*I unique Level 2 ID's in total)
#E.g. I=2 clinics, J=3 physicians per clinic: total number of unique physician ID's goes from 1-6
#n = no. of individuals, Level 1
#Total observations in dataset: N = I x J x n

#Takes in seed, number of Level 3 units, number of Level 2 units, number of Level 1 units, as well as:
#variances for random intercepts incorporated at Level 3 and Level 2, respectively; 
#variance for error included in treatment model; 
#and fixed effect coefficients (alpha.0-alpha.2)

data.gen <- function(seed, I, J, n, tau3, tau2, sigma2, alpha.0, alpha.1, alpha.2){

  #Set seed for reproducibility
  set.seed(seed)
  
  #Set total no. of observations
  N <- I*J*n
  
  #Generate Level 3 Column
  L3 <- rep(1:I, each=N/I)
  
  #Generate Level 2 column
  L2 <- rep(1:(J*I), each=n)

  #Initialize empty dataframe; one row for each individual
  trt.df <- matrix(, nrow=N, ncol=14)

  index <- 1

  #Iterate through each hierarchy (i.e. each Level 3 group)
  for(i in 1:I){

    #Generate random intercept for Level 3
    beta.0i <- sqrt(tau3)*rnorm(1, 0, 1)

    #Generate Stage 1 cluster-level covariate for Level 3
    V.i.1 <- rnorm(1, 1, 1)

    #Generate Stage 2 cluster-level covariate for Level 3
    V.i.2 <- rnorm(1, 1+0.5*V.i.1,  1)

    #Iterate through each Level 2 group
    for(j in 1:J){

      #Generate random intercept for Level 2
      beta.0j <- sqrt(tau2)*rnorm(1, 0, 1)

      #Iterate through each unit in Level 2
      for(k in 1:n){

        #Generate an individual level covariate for Stage 1
        X.ijk.1 <- rnorm(1, 2, 1)

        #Generate an individual level covariate for Stage 2
        X.ijk.2 <- rnorm(1, X.ijk.1, 1)

        #Generate error for treatment model
        eps.ijk.1 <- sqrt(sigma2)*rnorm(1, 0, 1)
        eps.ijk.2 <- sqrt(sigma2)*rnorm(1, 0, 1)

        #Treatment model with random intercepts: this is the TRUE UNDERLYING TREATMENT MODEL
        alpha.3 <- 0.5*alpha.2

        logit.p.1 <- alpha.1*X.ijk.1 + alpha.2*V.i.1 + alpha.3*(sin(X.ijk.1)) + alpha.0 + beta.0i + beta.0j + eps.ijk.1

        logit.p.2 <- alpha.1*X.ijk.2 + alpha.2*V.i.2 + alpha.3*(sin(X.ijk.2)) + alpha.0 + beta.0i + beta.0j + eps.ijk.2

        #Treatment model without random intercepts (hence "fixed")
        logit.p.fixed.1 <- alpha.1*X.ijk.1 + alpha.2*V.i.1 + alpha.3*(sin(X.ijk.1)) + alpha.0 + eps.ijk.1

        logit.p.fixed.2 <- alpha.1*X.ijk.2 + alpha.2*V.i.2 + alpha.3*(sin(X.ijk.2)) + alpha.0 + eps.ijk.2

        #Probability of treatment (trt) is based on true underlying treatment model
        prob.trt.1 <- exp(logit.p.1)/(exp(logit.p.1)+1)
        A.ijk.1 <- rbinom(1, 1, prob.trt.1)

        prob.trt.2 <- exp(logit.p.2)/(exp(logit.p.2)+1)
        A.ijk.2 <- rbinom(1, 1, prob.trt.2)

        #Indicators for whether an individual receives treatment or not
        #Not sure if we need this later on
        trt.indicator.1 <- ifelse(A.ijk.1 == 1, 1, -1)
        trt.indicator.2 <- ifelse(A.ijk.2 == 1, 1, -1)

        #Store data
        trt.df[index, 1:14] <- c(beta.0i, V.i.1, beta.0j, X.ijk.1, trt.indicator.1, logit.p.fixed.1, logit.p.1, A.ijk.1, V.i.2, X.ijk.2, trt.indicator.2, logit.p.fixed.2, logit.p.2, A.ijk.2)
        
        index <- index + 1
      }
    }
  }

  trt.data <- cbind(L3, L2, trt.df)
  trt.data <- data.frame(trt.data)
  trt.data$ID <- seq(1, N, 1)

  names(trt.data) <- c("L3", "L2", "beta.0i", "V.i.1", "beta.0j", "X.ijk.1", "trt.indicator.1", "logit.p.fixed.1", "logit.p.1", "A.ijk.1", "V.i.2", "X.ijk.2", "trt.indicator.2", "logit.p.fixed.2", "logit.p.2", "A.ijk.2", "ID")

  #MIS-SPECIFIED TREATMENT MODEL IGNORES NON-LINEAR TERM IN X.ijk

  #Fixed-effects only model, stage 1
  mod.1 <- glm(A.ijk.1~X.ijk.1 + I(sin(X.ijk.1)) + V.i.1, family = 'binomial', data = trt.data)
  
  #Fixed- and random-effects model, stage 1
  mod.2 <- glmer(A.ijk.1~X.ijk.1 + I(sin(X.ijk.1)) + V.i.1 + (1 | L3/L2), family = 'binomial', data = trt.data)

  #Fixed-effects only model, stage 2
  mod.3 <- glm(A.ijk.2~X.ijk.2 + I(sin(X.ijk.2)) + V.i.2, family = 'binomial', data = trt.data)
  
  #Fixed- and random-effects model, stage 2
  mod.4 <- glmer(A.ijk.2~X.ijk.2 + I(sin(X.ijk.2)) + V.i.2 + (1 | L3/L2), family = 'binomial', data = trt.data)
 
  #Store the estimated standard deviations of the random effects
  #Assuming independence of the random effects, we can sum the standard deviations for the standard deviation of 
  #the normal distribution of the random effects
  rand.sd.1 <- sum(as.data.frame(VarCorr(mod.2))[1:2,5])

  rand.sd.2 <- sum(as.data.frame(VarCorr(mod.4))[1:2,5])

  #IMPORTANT: now, logit.p.fixed and logit.p are different from what we defined earlier
  #logit.p.fixed.1 and logit.p.fixed.2 are the ESTIMATED logits of probability of treatment, integrating out the random effects of mod.2 and mod.4, respectively
  #Refer to stackexchange for why we set re.form=NA
  #Integrating out the random effects seems to be consistent with how the joint propensity score with random effects
  #is defined in the Perez-Heydrich et al. (2014) paper
  #logit.p.1 and logit.p.2 are the ESTIMATED logits of probability of treatment using mod.1 and mod.3, respectively
  trt.data$logit.p.fixed.1 <- predict(mod.2, trt.data, type='response', re.form=NULL)
  
  trt.data$logit.p.1 <- predict(mod.1, trt.data, type='response')

  trt.data$logit.p.fixed.2 <- predict(mod.4, trt.data, type='response', re.form=NULL)
  
  trt.data$logit.p.2 <- predict(mod.3, trt.data, type='response')
  
  #Return a dataframe with N x 17 columns 
  return(list(trt.data=trt.data, rand.sd.1=rand.sd.1, rand.sd.2=rand.sd.2))
}

#**********EGO AND ALTER SELECTION********** 

#Takes seed, generated treatment data, number of Level 3 groups, number of Level 2 groups,
#number of units in each Level 2 group, and number of egos (sam.size) we want to select

egos.alters <- function(seed, trt.data, I, J, n, sam.size){

  set.seed(seed)

  mydata.egos <- list()
  mydata.alters <- list()

  #Subset data based on Level 3
  for(i in 1:I){
    
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

#**********OUTCOME GENERATION**********

#Psi.0-psi.3 are the blip parameters

egos.outcomes <- function(trt.data, I, J, n, psi.0, psi.1, psi.2, psi.3){
  #Rename columns 17 and 35 so that we're distinguishing which ID's correspond to
  #egos and which ID's correspond to alters
  names(trt.data)[17] <- "ID.Ego"
  names(trt.data)[35] <- "ID.Alter"
  names(trt.data)[36] <- "alter"

  #alters.joint() was initially a product of alters' propensities, later becomes total number of treated
  #alters (accounting for both within- and between-group status)
  #For generating treatment decision rules and outcomes, we want to keep 
  #track of total number of within-group and between-group neighbours an ego has.
  #alters.cov.1 and alters.cov.2 will be summaries of the alters' covariates (e.g. the mean of alters' covariates)
  #alterswithin.cov.1 and alterswithin.cov.2 will be summaries of the within-group alters' covariates (e.g. the mean of the within-group alters' covariates)
  alters.joint.1 <- c()
  total.w.1 <- c()
  total.b.1 <- c()
  alter.cov.1 <- c()
  alterswithin.cov.1 <- c()

  alters.joint.2 <- c()
  total.w.2 <- c()
  total.b.2 <- c()
  alter.cov.2 <- c()
  alterswithin.cov.2 <- c()

  #Iterate through egos, first through Level 3 and then through Level 2 in Level 3
  for(i in 1:I){

    df.L3 <- trt.data[trt.data$L3 == i,]
    ego.w.1 <- c()
    ego.b.1 <- c()
    mean.alter.1 <- c()
    mean.alterswithin.1 <- c()

    ego.w.2 <- c()
    ego.b.2 <- c()
    mean.alter.2 <- c()
    mean.alterswithin.2 <- c()

    #This will be used for alters.joint().
    alters.prop.prod.1 <- c()

    alters.prop.prod.2 <- c()

    #We need to search through each unique ego ID in our dataframe (since
    #there are multiple rows for an ego, each linked with an alters'information).
    unique.egos <- unique(df.L3$ID.Ego)
    for(u in unique.egos){

      #Extract dataframe for that ego
      df.ego <- df.L3[df.L3$ID.Ego==u,]

      #Find its within-group neighbours
      df.ego.w <- df.ego[df.ego$alter.group == "W",]

      sum.trt.w.1 <- sum(df.ego.w[,28])

      sum.trt.w.2 <- sum(df.ego.w[,34])

      mean.alterswithin.1 <- c(mean.alterswithin.1, rep(mean(df.ego.w[,24]), nrow(df.ego)))
      mean.alterswithin.2 <- c(mean.alterswithin.2, rep(mean(df.ego.w[,30]), nrow(df.ego)))

      #Find its between-group neighbours
      df.ego.b <- df.ego[df.ego$alter.group == "B",]

      sum.trt.b.1 <- sum(df.ego.b[,28])

      sum.trt.b.2 <- sum(df.ego.b[,34])

      ego.w.1 <- c(ego.w.1, rep(sum.trt.w.1, nrow(df.ego)))

      ego.b.1 <- c(ego.b.1, rep(sum.trt.b.1, nrow(df.ego)))

      ego.w.2 <- c(ego.w.2, rep(sum.trt.w.2, nrow(df.ego)))

      ego.b.2 <- c(ego.b.2, rep(sum.trt.b.2, nrow(df.ego)))

      #Take the product of the ego's alters' propensities
      alters.prop.prod.1 <- c(alters.prop.prod.1, rep(prod(df.ego[,27]), nrow(df.ego)))
      alters.prop.prod.2 <- c(alters.prop.prod.2, rep(prod(df.ego[,33]), nrow(df.ego)))

      mean.alter.1 <- c(mean.alter.1, rep(mean(df.ego[,24]), nrow(df.ego)))
      mean.alter.2 <- c(mean.alter.2, rep(mean(df.ego[,30]), nrow(df.ego)))

    }

  alters.joint.1 <- c(alters.joint.1, alters.prop.prod.1)
  alters.joint.2 <- c(alters.joint.2, alters.prop.prod.2)
  total.w.1 <- c(total.w.1, ego.w.1)
  total.b.1 <- c(total.b.1, ego.b.1)
  total.w.2 <- c(total.w.2, ego.w.2)
  total.b.2 <- c(total.b.2, ego.b.2)
  alter.cov.1 <- c(alter.cov.1, mean.alter.1)
  alter.cov.2 <- c(alter.cov.2, mean.alter.2)

  alterswithin.cov.1 <- c(alterswithin.cov.1, mean.alterswithin.1)
  alterswithin.cov.2 <- c(alterswithin.cov.2, mean.alterswithin.2)
  }
  
  trt.data$alters.joint.1 <- alters.joint.1
  trt.data$alters.joint.2 <- alters.joint.2
  trt.data$total.w.1 <- total.w.1
  trt.data$total.b.1 <- total.b.1
  trt.data$total.w.2 <- total.w.2
  trt.data$total.b.2 <- total.b.2
  #This is where alters.joint.1 and alters.joint.2 turn into the total number of treated neighbours
  trt.data$alters.joint.1 <- total.w.1 + total.b.1
  trt.data$alters.joint.2 <- total.w.2 + total.b.2
  trt.data$alter.cov.1 <- alter.cov.1
  trt.data$alter.cov.2 <- alter.cov.2
  trt.data$alterswithin.cov.1 <- alterswithin.cov.1
  trt.data$alterswithin.cov.2 <- alterswithin.cov.2

  colnames(trt.data)[which(names(trt.data) == "alters.joint.1")] <- "total.t.1"
  colnames(trt.data)[which(names(trt.data) == "alters.joint.2")] <- "total.t.2"

  Y.ijk <- c()
  Y.err.ijk <- c()
  AOpt.ijk.1 <- c()
  AOpt.ijk.2 <- c()

  #Iterate through Level 3 groups
  for(i in 1:I){

    #Extract ego-specific information needed
    df.L3 <- trt.data[trt.data$L3 == i,]
    Y <- c()
    AOpt.1 <- c()
    AOpt.2 <- c()
    Y.err <- c()
    #Iterate through unique ego ID's
    unique.egos <- unique(df.L3$ID.Ego)
    for(u in unique.egos){
      #Extract a given ego's information (all of their rows, this is fine
      #since their cluster- and patient-specific information are the same across rows).
      df.ego <- df.L3[df.L3$ID.Ego==u,]

      #Treatment-free model: comprised of individual and cluster-level covariates
      #Interference terms account for both within- and between-groups of neighbours (i.e. total.t)
      #Interference term is interacting with cluster-level covariate V.i
      trt.free <- log(abs(df.ego$V.i.1)) + log(abs(df.ego$V.i.2)) + exp(df.ego$X.ijk.1) + exp(df.ego$X.ijk.2) + abs(df.ego$alter.cov.1) + abs(df.ego$alter.cov.2) + df.ego$total.t.1 + df.ego$total.t.1*df.ego$V.i.1 + df.ego$total.t.2 + df.ego$total.t.2*df.ego$V.i.2

      #Optimal treatment rules at stages 1 and 2 based on blip components
      A.Opt.1 <- as.numeric((df.ego$A.ijk.1*(psi.0 + psi.1*df.ego$X.ijk.1 + psi.2*df.ego$total.t.1 + psi.3*df.ego$V.i.1*df.ego$total.t.1))>0)
      A.Opt.2 <- as.numeric((df.ego$A.ijk.2*(psi.0 + psi.1*df.ego$X.ijk.2 + psi.2*df.ego$total.t.2 + psi.3*df.ego$V.i.2*df.ego$total.t.2))>0)

      mu1 <- (A.Opt.1-df.ego$A.ijk.1)*(psi.0 + psi.1*df.ego$X.ijk.1 + psi.2*df.ego$total.t.1 + psi.3*df.ego$V.i.1*df.ego$total.t.1)

      mu2 <- (A.Opt.2-df.ego$A.ijk.2)*(psi.0 + psi.1*df.ego$X.ijk.2 + psi.2*df.ego$total.t.2 + psi.3*df.ego$V.i.2*df.ego$total.t.2)

      blip1 <- (df.ego$A.ijk.1)*(psi.0 + psi.1*df.ego$X.ijk.1 + psi.2*df.ego$total.t.1 + psi.3*df.ego$V.i.1*df.ego$total.t.1)

      blip2 <- (df.ego$A.ijk.2)*(psi.0 + psi.1*df.ego$X.ijk.2 + psi.2*df.ego$total.t.2 + psi.3*df.ego$V.i.2*df.ego$total.t.2)

      #Overall outcome model
      mu.ijk <- trt.free + blip1 + blip2

      #Error component for outcome model
      y.err <- rep(rnorm(1, 0, 1), length(mu.ijk))
      #Y.out <- rnorm(length(mu.ijk), mu.ijk, 1)
      Y.out <- rep(rnorm(1, mu.ijk[1], 1), length(mu.ijk))

      Y <- c(Y, Y.out)
      AOpt.1 <- c(AOpt.1, A.Opt.1)
      AOpt.2 <- c(AOpt.2, A.Opt.2)
      Y.err <- c(Y.err, y.err)
    }

    Y.ijk <- c(Y.ijk, Y)
    AOpt.ijk.1 <- c(AOpt.ijk.1, AOpt.1)
    AOpt.ijk.2 <- c(AOpt.ijk.2, AOpt.2)
    Y.err.ijk <- c(Y.err.ijk, Y.err)
  }

  trt.data$Y.ijk <- Y.ijk
  trt.data$AOpt.ijk.1 <- AOpt.ijk.1
  trt.data$AOpt.ijk.2 <- AOpt.ijk.2
  trt.data$Y.err.ijk <- Y.err.ijk

  return(trt.data)
}

#**********WEIGHT CONSTRUCTION**********

#The group.propensity function will compute different weights 
#used in our different methods

group.propensity <- function(trt.data, I, J, n, tau3, tau2, rand.sd.1, rand.sd.2){

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

    unique.egos <- unique(df.L3$ID.Ego)
    
    #Iterate through each ego within Level 3 (possibly
    #from different Level 2 groups)
    for(u in unique.egos){

      #Get dataframe corresponding to ego
      df.ego <- df.L3[df.L3$ID.Ego==u,]

      #Get the ego's probability of treatment from stage 1 mixed model (it's called logit but is actually the expit)
      #and indicator for whether or not they received treatment "1"
      ego.logit.1 <- df.ego$logit.p.fixed.1[1]
      ego.indicator.1 <- df.ego$trt.indicator.1[1]
      ego.A.1 <- ifelse(ego.indicator.1==1, 1, 0)

      #Get vector of alters' (within- and between-group) probability of treatment from stage 1 mixed model
      alter.logit.1 <- df.ego$logit.p.fixed.1.alter
      alter.indicator.1 <- df.ego$trt.indicator.1.alter
      alter.A.1 <- ifelse(alter.indicator.1==1, 1, 0)

      #Get the ego's probability of treatment from stage 2 mixed model (it's called logit but is actually the expit)
      #and indicator for whether or not they received treatment "1"
      ego.logit.2 <- df.ego$logit.p.fixed.2[1]
      ego.indicator.2 <- df.ego$trt.indicator.2[1]
      ego.A.2 <- ifelse(ego.indicator.2==1, 1, 0)

      #Get vector of alters' (within- and between-group) probability of treatment from stage 2 mixed model
      alter.logit.2 <- df.ego$logit.p.fixed.2.alter
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
        ego.w.logit.1 <- df.ego$logit.p.fixed.1[1]
        ego.w.indicator.1 <- df.ego$trt.indicator.1[1]
        ego.w.A.1 <- ifelse(ego.w.indicator.1==1, 1, 0)

        ego.w.logit.2 <- df.ego$logit.p.fixed.2[1]
        ego.w.indicator.2 <- df.ego$trt.indicator.2[1]
        ego.w.A.2 <- ifelse(ego.w.indicator.2==1, 1, 0)

        w.within.L2.1 <- c(w.within.L2.1, rep(abs(ego.w.A.1-ego.w.logit.1), nrow(df.ego)))
        w.network.L2.1 <- c(w.network.L2.1, rep(abs(ego.w.A.1-ego.w.logit.1), nrow(df.ego)))

        w.within.L2.2 <- c(w.within.L2.2, rep(abs(ego.w.A.2-ego.w.logit.2), nrow(df.ego)))
        w.network.L2.2 <- c(w.network.L2.2, rep(abs(ego.w.A.2-ego.w.logit.2), nrow(df.ego)))

      }

      else{
        ego.w.logit.1 <- df.ego$logit.p.fixed.1[1]
        ego.w.indicator.1 <- df.ego$trt.indicator.1[1]
        ego.w.A.1 <- ifelse(ego.w.indicator.1==1, 1, 0)

        alter.w.logit.1 <- df.ego.w$logit.p.fixed.1.alter
        alter.w.indicator.1 <- df.ego.w$trt.indicator.1.alter
        alter.w.A.1 <- ifelse(alter.w.indicator.1==1, 1, 0)

        ego.w.logit.2 <- df.ego$logit.p.fixed.2[1]
        ego.w.indicator.2 <- df.ego$trt.indicator.2[1]
        ego.w.A.2 <- ifelse(ego.w.indicator.2==1, 1, 0)

        alter.w.logit.2 <- df.ego.w$logit.p.fixed.2.alter
        alter.w.indicator.2 <- df.ego.w$trt.indicator.2.alter
        alter.w.A.2 <- ifelse(alter.w.indicator.2==1, 1, 0)

        w.within.L2.1 <- c(w.within.L2.1, rep(abs(ego.w.A.1-ego.w.logit.1)*prod(abs(alter.w.A.1-alter.w.logit.1)), nrow(df.ego)))
        w.network.L2.1 <- c(w.network.L2.1, rep(abs(ego.w.A.1-ego.w.logit.1)*prod(abs(alter.w.A.1-alter.w.logit.1)), nrow(df.ego)))

        w.within.L2.2 <- c(w.within.L2.2, rep(abs(ego.w.A.2-ego.w.logit.2)*prod(abs(alter.w.A.2-alter.w.logit.2)), nrow(df.ego)))
        w.network.L2.2 <- c(w.network.L2.2, rep(abs(ego.w.A.2-ego.w.logit.2)*prod(abs(alter.w.A.2-alter.w.logit.2)), nrow(df.ego)))

      }
      
      #Weights for standard individual absolute value weight
      #from mixed-effects model
      ego.s.logit.1 <- df.ego$logit.p.fixed.1[1]
      ego.s.indicator.1 <- df.ego$trt.indicator.1[1]
      ego.s.A.1 <- ifelse(ego.s.indicator.1==1, 1, 0)

      ipw.s.1 <- abs(ego.s.A.1 - ego.s.logit.1)
      w.standard.L2.1 <- c(w.standard.L2.1, rep(ipw.s.1, nrow(df.ego)))

      ego.s.logit.2 <- df.ego$logit.p.fixed.2[1]
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

#**********GENERATE DATASET FOR EGOS**********

#Initializes hierarchical network structure, selects egos and alters,
#generates outcomes, and calculates different weights
main.fun <- function(I, J, n, sam.size, tau3, tau2, sigma2, alpha.0, alpha.1, alpha.2, psi.0, psi.1, psi.2, psi.3, seed){
  data.gen.list <- data.gen(seed, I, J, n, tau3, tau2, sigma2, alpha.0, alpha.1, alpha.2)
  trt.data <- data.gen.list[[1]]
  rand.sd.1 <- data.gen.list[[2]]
  rand.sd.2 <- data.gen.list[[3]]
  trt.data <- egos.alters(seed, trt.data, I, J, n, sam.size)
  trt.data <- egos.outcomes(trt.data, I, J, n, psi.0, psi.1, psi.2, psi.3)
  trt.data <- group.propensity(trt.data, I, J, n, tau3, tau2, rand.sd.1, rand.sd.2)
  
  ego.df <- trt.data

  num.reps <- as.numeric(table(ego.df$ID.Ego))

  ego.df <- ego.df[order(ego.df$ID.Ego),]

  ego.df$num.alters <- rep(num.reps, num.reps) 

  ego.df <- ego.df[!duplicated(ego.df$ID.Ego), ]
  
  ego.df <- subset(ego.df, select=c(ID.Ego, L3:num.alters))

  ego.df <- ego.df[order(ego.df$num.alters),]

  w.within.mean.2 <- c()

  w.network.mean.2 <- c()

  w.L3.mean.2 <- c()

  w.within.mean.1 <- c()

  w.network.mean.1 <- c()

  w.L3.mean.1 <- c()

  
  for(t in unique(ego.df$num.alters)){

      df.t.2 <- ego.df[which(ego.df$num.alters==t),]
      w.within.mean.2 <- c(w.within.mean.2, rep(mean(df.t.2$w.within.2), nrow(df.t.2)))
      #print(mean(df.t$w.within))
      w.network.mean.2 <- c(w.network.mean.2, rep(mean(df.t.2$w.network.2), nrow(df.t.2)))
      #print(mean(df.t$w.network))
      w.L3.mean.2 <- c(w.L3.mean.2, rep(mean(df.t.2$w.L3.2), nrow(df.t.2)))
      #print(mean(df.t$w.L3))

      #df.t.1 <- ego.df[which(ego.df$num.alters==t),]
      w.within.mean.1 <- c(w.within.mean.1, rep(mean(df.t.2$w.within.1), nrow(df.t.2)))
      #print(mean(df.t$w.within))
      w.network.mean.1 <- c(w.network.mean.1, rep(mean(df.t.2$w.network.1), nrow(df.t.2)))
      #print(mean(df.t$w.network))
      w.L3.mean.1 <- c(w.L3.mean.1, rep(mean(df.t.2$w.L3.1), nrow(df.t.2)))
      #print(mean(df.t$w.L3))
  }
  

  ego.df$w.within.mean.1 <- w.within.mean.1
  ego.df$w.network.mean.1 <- w.network.mean.1
  ego.df$w.L3.mean.1 <- w.L3.mean.1

  ego.df$w.within.mean.2 <- w.within.mean.2
  ego.df$w.network.mean.2 <- w.network.mean.2
  ego.df$w.L3.mean.2 <- w.L3.mean.2
  
  #print(w.within.mean)
  #print(w.network.mean)
  #print(w.L3.mean)

  #w.within.mean.2 <- replace(w.within.mean.2, w.within.mean.2 == 0, 0.01)

  #w.network.mean.2 <- replace(w.network.mean.2, w.network.mean.2 == 0, 0.01)

  #w.L3.mean.2 <- replace(w.L3.mean.2, w.L3.mean.2 == 0, 0.01)


  #as.numeric(table(ego.df$total.t))
  #w.w.mean.2 <- rep(w.within.mean.2, times=as.numeric(table(ego.df$num.alters)))

  #w.n.mean.2 <- rep(w.network.mean.2, times=as.numeric(table(ego.df$num.alters)))

  #w.L.mean.2 <- rep(w.L3.mean.2, times=as.numeric(table(ego.df$num.alters)))

  ego.df$w.within.resc.2 <- ego.df$w.within.2/w.within.mean.2

  ego.df$w.network.resc.2 <- ego.df$w.network.2/w.network.mean.2

  ego.df$w.L3.resc.2 <- ego.df$w.L3.2/w.L3.mean.2


  #***********************************************************************************
  #ego.df <- ego.df[order(ego.df$total.t.1),]
  
  #print(w.within.mean)
  #print(w.network.mean)
  #print(w.L3.mean)

  #w.within.mean.1 <- replace(w.within.mean.1, w.within.mean.1 == 0, 0.01)

  #w.network.mean.1 <- replace(w.network.mean.1, w.network.mean.1 == 0, 0.01)

  #w.L3.mean.1 <- replace(w.L3.mean.1, w.L3.mean.1 == 0, 0.01)


  #as.numeric(table(ego.df$total.t))
  #w.w.mean.1 <- rep(w.within.mean.1, times=as.numeric(table(ego.df$num.alters)))

  #w.n.mean.1 <- rep(w.network.mean.1, times=as.numeric(table(ego.df$num.alters)))

  #w.L.mean.1 <- rep(w.L3.mean.1, times=as.numeric(table(ego.df$num.alters)))

  ego.df$w.within.resc.1 <- ego.df$w.within.1/w.within.mean.1

  ego.df$w.network.resc.1 <- ego.df$w.network.1/w.network.mean.1

  ego.df$w.L3.resc.1 <- ego.df$w.L3.1/w.L3.mean.1
  
  return(ego.df)
}

#Test out a few different datasets
ego.df1 <- main.fun(10, 4, 5, 2, 1, 1, 1, -1, 1, 0.5, 1, 1, -1, -1.5, 1000)
#sum(abs(ego.df1$AOpt.ijk-ego.df1$A.ijk))/nrow(ego.df1)

ego.df2 <- main.fun(10, 5, 5, 2, 1, 1, 1, -1, 1, 0.5, 1, 1, -1, -1.5, 1000)
#sum(abs(ego.df2$AOpt.ijk-ego.df2$A.ijk))/nrow(ego.df2)

ego.df3 <- main.fun(10, 6, 5, 2, 1, 1, 1, -1, 1, 0.5, 1, 1, -1, -1.5, 1000)
#sum(abs(ego.df3$AOpt.ijk-ego.df3$A.ijk))/nrow(ego.df3)

#**********DWOLS IMPLEMENTATION WHERE TREATMENT MODEL IS MIS-SPECIFIED WITH NO CLUSTER-LEVEL COVARIATE********** 

#ME4 AT BOTH STAGES

dwols.wthn <- function(data.sim, data.test, I, J, n) {
  N <- I*J*n

  mod2 <- lm(Y.ijk ~ V.i.1 + V.i.2 + X.ijk.1 + X.ijk.2 + I(abs(alterswithin.cov.1)) + I(abs(alterswithin.cov.2)) + total.w.1 + total.w.1:V.i.1 + total.w.2 + total.w.2:V.i.2 + A.ijk.1 + A.ijk.1:X.ijk.1 + A.ijk.1:total.w.1 + A.ijk.1:total.w.1:V.i.1 + A.ijk.2 + A.ijk.2:X.ijk.2 + A.ijk.2:total.w.2 + A.ijk.2:total.w.2:V.i.2, weights = w.within.2, data = data.sim)
  
  #Obtain extra blip parameter estimates from outcome model
  #(only coefficients for treatment terms)
  psi.0.stage2 <- coef(mod2)[11]
  psi.1.stage2 <- coef(mod2)[16]
  psi.2.stage2 <- coef(mod2)[17]
  psi.3.stage2 <- coef(mod2)[19]

  xsi.0.stage2 <- coef(mod2)[9]
  xsi.1.stage2 <- coef(mod2)[13]
  #Compute optimal treatment in testing set using blip parameter estimates
  Aopt.dwols.2 <- as.numeric((psi.0.stage2 + psi.1.stage2*data.sim$X.ijk.2 + psi.2.stage2*data.sim$total.w.2 + psi.3.stage2*data.sim$V.i.2*data.sim$total.w.2) > 0)
  
  #Optimal outcome
  data.sim$Y.ijk.2 <- data.sim$Y.ijk + (Aopt.dwols.2 - data.sim$A.ijk.2)*(psi.0.stage2 + psi.1.stage2*data.sim$X.ijk.2 + psi.2.stage2*data.sim$total.w.2 + psi.3.stage2*data.sim$V.i.2*data.sim$total.w.2)
  
  mod1 <- lm(Y.ijk.2 ~ V.i.1 + X.ijk.1 + I(abs(alterswithin.cov.1)) + total.w.1 + total.w.1:V.i.1 + A.ijk.1 + A.ijk.1:X.ijk.1 + A.ijk.1:total.w.1 + A.ijk.1:total.w.1:V.i.1, weights = w.within.1, data = data.sim)

  psi.0.stage1 <- coef(mod1)[6]
  psi.1.stage1 <- coef(mod1)[8]
  psi.2.stage1 <- coef(mod1)[9]
  psi.3.stage1 <- coef(mod1)[10]

  xsi.0.stage1 <- coef(mod1)[5]
  xsi.1.stage1 <- coef(mod1)[7]
  #E[Y|A = Aopt.dwols]
  #Optimal treatment based on what we observe in given dataset

  Aopt.dwols.1 <- as.numeric((psi.0.stage1 + psi.1.stage1*data.sim$X.ijk.1 + psi.2.stage1*data.sim$total.w.1 + psi.3.stage1*data.sim$V.i.1*data.sim$total.w.1) > 0)


  #Optimal treatment in testing set
  Aopt.dwols.2.test <- as.numeric((psi.0.stage2 + psi.1.stage2*data.test$X.ijk.2 + psi.2.stage2*data.test$total.w.2 + psi.3.stage2*data.test$V.i.2*data.test$total.w.2) > 0)
  Aopt.dwols.1.test <- as.numeric((psi.0.stage1 + psi.1.stage1*data.test$X.ijk.1 + psi.2.stage1*data.test$total.w.1 + psi.3.stage1*data.test$V.i.1*data.test$total.w.1) > 0)

  Yopt.test <- data.test$Y.ijk + (Aopt.dwols.2.test-data.test$A.ijk.2)*(psi.0.stage2 + psi.1.stage2*data.test$X.ijk.2 + psi.2.stage2*data.test$total.w.2 + psi.3.stage2*data.test$V.i.2*data.test$total.w.2) + (Aopt.dwols.1.test-data.test$A.ijk.1)*(psi.0.stage1 + psi.1.stage1*data.test$X.ijk.1 + psi.2.stage1*data.test$total.w.1 + psi.3.stage1*data.test$V.i.1*data.test$total.w.1)
  #Aopt.dwols.obs <- as.numeric((psi.0.hat + psi.1.hat*data.sim$X.ijk) > 0)
  
  #Mean of outcomes where treatment we observe is the same as the dwols optimal
  #treatment for our given dataset
  #Y.agree <- mean(data.sim$Y.ijk[data.sim$A.ijk == Aopt.dwols.obs])
  
  return(data.frame(Aopt.dwols.1.test, Aopt.dwols.2.test, Yopt.test, psi.0.stage1, psi.1.stage1, psi.2.stage1, psi.3.stage1, psi.0.stage2, psi.1.stage2, psi.2.stage2, psi.3.stage2, xsi.0.stage1, xsi.1.stage1, xsi.0.stage2, xsi.1.stage2))
}

dwols.standard <- function(data.sim, data.test, I, J, n) {
  N <- I*J*n

  mod2 <- lm(Y.ijk ~ V.i.1 + V.i.2 + X.ijk.1 + X.ijk.2 + I(abs(alter.cov.1)) + I(abs(alter.cov.2)) + total.t.1 + total.t.1:V.i.1 + total.t.2 + total.t.2:V.i.2 + A.ijk.1 + A.ijk.1:X.ijk.1 + A.ijk.1:total.t.1 + A.ijk.1:total.t.1:V.i.1 + A.ijk.2 + A.ijk.2:X.ijk.2 + A.ijk.2:total.t.2 + A.ijk.2:total.t.2:V.i.2, weights = w.standard.2, data = data.sim)
  
  #Obtain extra blip parameter estimates from outcome model
  #(only coefficients for treatment terms)
  psi.0.stage2 <- coef(mod2)[11]
  psi.1.stage2 <- coef(mod2)[16]
  psi.2.stage2 <- coef(mod2)[17]
  psi.3.stage2 <- coef(mod2)[19]

  xsi.0.stage2 <- coef(mod2)[9]
  xsi.1.stage2 <- coef(mod2)[13]
  #Compute optimal treatment in testing set using blip parameter estimates
  Aopt.dwols.2 <- as.numeric((psi.0.stage2 + psi.1.stage2*data.sim$X.ijk.2 + psi.2.stage2*data.sim$total.t.2 + psi.3.stage2*data.sim$V.i.2*data.sim$total.t.2) > 0)
  
  #Optimal outcome
  data.sim$Y.ijk.2 <- data.sim$Y.ijk + (Aopt.dwols.2 - data.sim$A.ijk.2)*(psi.0.stage2 + psi.1.stage2*data.sim$X.ijk.2 + psi.2.stage2*data.sim$total.t.2 + psi.3.stage2*data.sim$V.i.2*data.sim$total.t.2)
  
  mod1 <- lm(Y.ijk.2 ~ V.i.1 + X.ijk.1 + I(abs(alter.cov.1)) + total.t.1 + total.t.1:V.i.1 + A.ijk.1 + A.ijk.1:X.ijk.1 + A.ijk.1:total.t.1 + A.ijk.1:total.t.1:V.i.1, weights = w.standard.1, data = data.sim)

  psi.0.stage1 <- coef(mod1)[6]
  psi.1.stage1 <- coef(mod1)[8]
  psi.2.stage1 <- coef(mod1)[9]
  psi.3.stage1 <- coef(mod1)[10]

  xsi.0.stage1 <- coef(mod1)[5]
  xsi.1.stage1 <- coef(mod1)[7]
  #E[Y|A = Aopt.dwols]
  #Optimal treatment based on what we observe in given dataset

  Aopt.dwols.1 <- as.numeric((psi.0.stage1 + psi.1.stage1*data.sim$X.ijk.1 + psi.2.stage1*data.sim$total.t.1 + psi.3.stage1*data.sim$V.i.1*data.sim$total.t.1) > 0)


  #Optimal treatment in testing set
  Aopt.dwols.2.test <- as.numeric((psi.0.stage2 + psi.1.stage2*data.test$X.ijk.2 + psi.2.stage2*data.test$total.t.2 + psi.3.stage2*data.test$V.i.2*data.test$total.t.2) > 0)
  Aopt.dwols.1.test <- as.numeric((psi.0.stage1 + psi.1.stage1*data.test$X.ijk.1 + psi.2.stage1*data.test$total.t.1 + psi.3.stage1*data.test$V.i.1*data.test$total.t.1) > 0)

  Yopt.test <- data.test$Y.ijk + (Aopt.dwols.2.test-data.test$A.ijk.2)*(psi.0.stage2 + psi.1.stage2*data.test$X.ijk.2 + psi.2.stage2*data.test$total.t.2 + psi.3.stage2*data.test$V.i.2*data.test$total.t.2) + (Aopt.dwols.1.test-data.test$A.ijk.1)*(psi.0.stage1 + psi.1.stage1*data.test$X.ijk.1 + psi.2.stage1*data.test$total.t.1 + psi.3.stage1*data.test$V.i.1*data.test$total.t.1)
  #Aopt.dwols.obs <- as.numeric((psi.0.hat + psi.1.hat*data.sim$X.ijk) > 0)
  
  #Mean of outcomes where treatment we observe is the same as the dwols optimal
  #treatment for our given dataset
  #Y.agree <- mean(data.sim$Y.ijk[data.sim$A.ijk == Aopt.dwols.obs])
  
  return(data.frame(Aopt.dwols.1.test, Aopt.dwols.2.test, Yopt.test, psi.0.stage1, psi.1.stage1, psi.2.stage1, psi.3.stage1, psi.0.stage2, psi.1.stage2, psi.2.stage2, psi.3.stage2, xsi.0.stage1, xsi.1.stage1, xsi.0.stage2, xsi.1.stage2))
}

dwols.both <- function(data.sim, data.test, I, J, n) {
  N <- I*J*n

  mod2 <- lm(Y.ijk ~ V.i.1 + V.i.2 + X.ijk.1 + X.ijk.2 + I(abs(alter.cov.1)) + I(abs(alter.cov.2)) + total.t.1 + total.t.1:V.i.1 + total.t.2 + total.t.2:V.i.2 + A.ijk.1 + A.ijk.1:X.ijk.1 + A.ijk.1:total.t.1 + A.ijk.1:total.t.1:V.i.1 + A.ijk.2 + A.ijk.2:X.ijk.2 + A.ijk.2:total.t.2 + A.ijk.2:total.t.2:V.i.2, weights = w.L3.2, data = data.sim)
  
  #Obtain extra blip parameter estimates from outcome model
  #(only coefficients for treatment terms)
  psi.0.stage2 <- coef(mod2)[11]
  psi.1.stage2 <- coef(mod2)[16]
  psi.2.stage2 <- coef(mod2)[17]
  psi.3.stage2 <- coef(mod2)[19]

  xsi.0.stage2 <- coef(mod2)[9]
  xsi.1.stage2 <- coef(mod2)[13]
  #Compute optimal treatment in testing set using blip parameter estimates
  Aopt.dwols.2 <- as.numeric((psi.0.stage2 + psi.1.stage2*data.sim$X.ijk.2 + psi.2.stage2*data.sim$total.t.2 + psi.3.stage2*data.sim$V.i.2*data.sim$total.t.2) > 0)
  
  #Optimal outcome
  data.sim$Y.ijk.2 <- data.sim$Y.ijk + (Aopt.dwols.2 - data.sim$A.ijk.2)*(psi.0.stage2 + psi.1.stage2*data.sim$X.ijk.2 + psi.2.stage2*data.sim$total.t.2 + psi.3.stage2*data.sim$V.i.2*data.sim$total.t.2)
  
  mod1 <- lm(Y.ijk.2 ~ V.i.1 + X.ijk.1 + I(abs(alter.cov.1)) + total.t.1 + total.t.1:V.i.1 + A.ijk.1 + A.ijk.1:X.ijk.1 + A.ijk.1:total.t.1 + A.ijk.1:total.t.1:V.i.1, weights = w.L3.1, data = data.sim)

  psi.0.stage1 <- coef(mod1)[6]
  psi.1.stage1 <- coef(mod1)[8]
  psi.2.stage1 <- coef(mod1)[9]
  psi.3.stage1 <- coef(mod1)[10]

  xsi.0.stage1 <- coef(mod1)[5]
  xsi.1.stage1 <- coef(mod1)[7]
  #E[Y|A = Aopt.dwols]
  #Optimal treatment based on what we observe in given dataset

  Aopt.dwols.1 <- as.numeric((psi.0.stage1 + psi.1.stage1*data.sim$X.ijk.1 + psi.2.stage1*data.sim$total.t.1 + psi.3.stage1*data.sim$V.i.1*data.sim$total.t.1) > 0)


  #Optimal treatment in testing set
  Aopt.dwols.2.test <- as.numeric((psi.0.stage2 + psi.1.stage2*data.test$X.ijk.2 + psi.2.stage2*data.test$total.t.2 + psi.3.stage2*data.test$V.i.2*data.test$total.t.2) > 0)
  Aopt.dwols.1.test <- as.numeric((psi.0.stage1 + psi.1.stage1*data.test$X.ijk.1 + psi.2.stage1*data.test$total.t.1 + psi.3.stage1*data.test$V.i.1*data.test$total.t.1) > 0)

  Yopt.test <- data.test$Y.ijk + (Aopt.dwols.2.test-data.test$A.ijk.2)*(psi.0.stage2 + psi.1.stage2*data.test$X.ijk.2 + psi.2.stage2*data.test$total.t.2 + psi.3.stage2*data.test$V.i.2*data.test$total.t.2) + (Aopt.dwols.1.test-data.test$A.ijk.1)*(psi.0.stage1 + psi.1.stage1*data.test$X.ijk.1 + psi.2.stage1*data.test$total.t.1 + psi.3.stage1*data.test$V.i.1*data.test$total.t.1)
  #Aopt.dwols.obs <- as.numeric((psi.0.hat + psi.1.hat*data.sim$X.ijk) > 0)
  
  #Mean of outcomes where treatment we observe is the same as the dwols optimal
  #treatment for our given dataset
  #Y.agree <- mean(data.sim$Y.ijk[data.sim$A.ijk == Aopt.dwols.obs])
  
  return(data.frame(Aopt.dwols.1.test, Aopt.dwols.2.test, Yopt.test, psi.0.stage1, psi.1.stage1, psi.2.stage1, psi.3.stage1, psi.0.stage2, psi.1.stage2, psi.2.stage2, psi.3.stage2, xsi.0.stage1, xsi.1.stage1, xsi.0.stage2, xsi.1.stage2))
}


#ME2 AT BOTH STAGES
dwols.both.resc <- function(data.sim, data.test, I, J, n) {
  N <- I*J*n

  mod2 <- lm(Y.ijk ~ V.i.1 + V.i.2 + X.ijk.1 + X.ijk.2 + I(abs(alter.cov.1)) + I(abs(alter.cov.2)) + total.t.1 + total.t.1:V.i.1 + total.t.2 + total.t.2:V.i.2 + A.ijk.1 + A.ijk.1:X.ijk.1 + A.ijk.1:total.t.1 + A.ijk.1:total.t.1:V.i.1 + A.ijk.2 + A.ijk.2:X.ijk.2 + A.ijk.2:total.t.2 + A.ijk.2:total.t.2:V.i.2, weights = w.L3.resc.2, data = data.sim)
  
  #Obtain extra blip parameter estimates from outcome model
  #(only coefficients for treatment terms)
  psi.0.stage2 <- coef(mod2)[11]
  psi.1.stage2 <- coef(mod2)[16]
  psi.2.stage2 <- coef(mod2)[17]
  psi.3.stage2 <- coef(mod2)[19]

  xsi.0.stage2 <- coef(mod2)[9]
  xsi.1.stage2 <- coef(mod2)[13]
  #Compute optimal treatment in testing set using blip parameter estimates
  Aopt.dwols.2 <- as.numeric((psi.0.stage2 + psi.1.stage2*data.sim$X.ijk.2 + psi.2.stage2*data.sim$total.t.2 + psi.3.stage2*data.sim$V.i.2*data.sim$total.t.2) > 0)
  
  #Optimal outcome
  data.sim$Y.ijk.2 <- data.sim$Y.ijk + (Aopt.dwols.2 - data.sim$A.ijk.2)*(psi.0.stage2 + psi.1.stage2*data.sim$X.ijk.2 + psi.2.stage2*data.sim$total.t.2 + psi.3.stage2*data.sim$V.i.2*data.sim$total.t.2)
  
  mod1 <- lm(Y.ijk.2 ~ V.i.1 + X.ijk.1 + I(abs(alter.cov.1)) + total.t.1 + total.t.1:V.i.1 + A.ijk.1 + A.ijk.1:X.ijk.1 + A.ijk.1:total.t.1 + A.ijk.1:total.t.1:V.i.1, weights = w.L3.resc.1, data = data.sim)

  psi.0.stage1 <- coef(mod1)[6]
  psi.1.stage1 <- coef(mod1)[8]
  psi.2.stage1 <- coef(mod1)[9]
  psi.3.stage1 <- coef(mod1)[10]

  xsi.0.stage1 <- coef(mod1)[5]
  xsi.1.stage1 <- coef(mod1)[7]
  #E[Y|A = Aopt.dwols]
  #Optimal treatment based on what we observe in given dataset

  Aopt.dwols.1 <- as.numeric((psi.0.stage1 + psi.1.stage1*data.sim$X.ijk.1 + psi.2.stage1*data.sim$total.t.1 + psi.3.stage1*data.sim$V.i.1*data.sim$total.t.1) > 0)


  #Optimal treatment in testing set
  Aopt.dwols.2.test <- as.numeric((psi.0.stage2 + psi.1.stage2*data.test$X.ijk.2 + psi.2.stage2*data.test$total.t.2 + psi.3.stage2*data.test$V.i.2*data.test$total.t.2) > 0)
  Aopt.dwols.1.test <- as.numeric((psi.0.stage1 + psi.1.stage1*data.test$X.ijk.1 + psi.2.stage1*data.test$total.t.1 + psi.3.stage1*data.test$V.i.1*data.test$total.t.1) > 0)

  Yopt.test <- data.test$Y.ijk + (Aopt.dwols.2.test-data.test$A.ijk.2)*(psi.0.stage2 + psi.1.stage2*data.test$X.ijk.2 + psi.2.stage2*data.test$total.t.2 + psi.3.stage2*data.test$V.i.2*data.test$total.t.2) + (Aopt.dwols.1.test-data.test$A.ijk.1)*(psi.0.stage1 + psi.1.stage1*data.test$X.ijk.1 + psi.2.stage1*data.test$total.t.1 + psi.3.stage1*data.test$V.i.1*data.test$total.t.1)
  #Aopt.dwols.obs <- as.numeric((psi.0.hat + psi.1.hat*data.sim$X.ijk) > 0)
  
  #Mean of outcomes where treatment we observe is the same as the dwols optimal
  #treatment for our given dataset
  #Y.agree <- mean(data.sim$Y.ijk[data.sim$A.ijk == Aopt.dwols.obs])
  
  return(data.frame(Aopt.dwols.1.test, Aopt.dwols.2.test, Yopt.test, psi.0.stage1, psi.1.stage1, psi.2.stage1, psi.3.stage1, psi.0.stage2, psi.1.stage2, psi.2.stage2, psi.3.stage2, xsi.0.stage1, xsi.1.stage1, xsi.0.stage2, xsi.1.stage2))
}

#******************************************************************************************************************************
#SIMULATIONS START HERE
Sims <- function(I, J, n, sam.size, tau3, tau2, sigma2, alpha.0, alpha.1, alpha.2, psi.0, psi.1, psi.2, psi.3, run, model) {
  Res.both <- data.frame(A1OPT = double(), A2OPT = double(), DWOLS.BOTH.STAGE1 = double(), DWOLS.BOTH.STAGE2 = double(),  Yopt = double(), Yopt.DWOLS.BOTH = double(), PSI.0.STAGE1.BOTH=double(), PSI.1.STAGE1.BOTH=double(),PSI.2.STAGE1.BOTH=double(),PSI.3.STAGE1.BOTH=double(),PSI.0.STAGE2.BOTH=double(), PSI.1.STAGE2.BOTH=double(),PSI.2.STAGE2.BOTH=double(),PSI.3.STAGE2.BOTH=double(), XSI.0.STAGE1.BOTH=double(), XSI.1.STAGE1.BOTH=double(), XSI.0.STAGE2.BOTH=double(), XSI.1.STAGE2.BOTH=double())
  Res.both.resc <- data.frame(A1OPT = double(), A2OPT = double(), DWOLS.BOTH.RESC.STAGE1 = double(), DWOLS.BOTH.RESC.STAGE2 = double(),  Yopt = double(), Yopt.DWOLS.BOTH.RESC = double(), PSI.0.STAGE1.BOTH.RESC=double(), PSI.1.STAGE1.BOTH.RESC=double(),PSI.2.STAGE1.BOTH.RESC=double(),PSI.3.STAGE1.BOTH.RESC=double(),PSI.0.STAGE2.BOTH.RESC=double(), PSI.1.STAGE2.BOTH.RESC=double(),PSI.2.STAGE2.BOTH.RESC=double(),PSI.3.STAGE2.BOTH.RESC=double(), XSI.0.STAGE1.BOTH.RESC=double(), XSI.1.STAGE1.BOTH.RESC=double(), XSI.0.STAGE2.BOTH.RESC=double(), XSI.1.STAGE2.BOTH.RESC=double())

  Res.wthn <- data.frame(A1OPT = double(), A2OPT = double(), DWOLS.WTHN.STAGE1 = double(), DWOLS.WTHN.STAGE2 = double(),  Yopt = double(), Yopt.DWOLS.WTHN = double(), PSI.0.STAGE1.WTHN=double(), PSI.1.STAGE1.WTHN=double(),PSI.2.STAGE1.WTHN=double(),PSI.3.STAGE1.WTHN=double(),PSI.0.STAGE2.WTHN=double(), PSI.1.STAGE2.WTHN=double(),PSI.2.STAGE2.WTHN=double(),PSI.3.STAGE2.WTHN=double(), XSI.0.STAGE1.WTHN=double(), XSI.1.STAGE1.WTHN=double(), XSI.0.STAGE2.WTHN=double(), XSI.1.STAGE2.WTHN=double())
  Res.standard <- data.frame(A1OPT = double(), A2OPT = double(), DWOLS.STANDARD.STAGE1 = double(), DWOLS.STANDARD.STAGE2 = double(),  Yopt = double(), Yopt.DWOLS.STANDARD = double(), PSI.0.STAGE1.STANDARD=double(), PSI.1.STAGE1.STANDARD=double(),PSI.2.STAGE1.STANDARD=double(),PSI.3.STAGE1.STANDARD=double(),PSI.0.STAGE2.STANDARD=double(), PSI.1.STAGE2.STANDARD=double(),PSI.2.STAGE2.STANDARD=double(),PSI.3.STAGE2.STANDARD=double(), XSI.0.STAGE1.STANDARD=double(), XSI.1.STAGE1.STANDARD=double(), XSI.0.STAGE2.STANDARD=double(), XSI.1.STAGE2.STANDARD=double())

  Min.Num.Alters <- c()
  Max.Num.Alters <- c()
  Mean.Num.Alters <- c()
  Median.Num.Alters <- c()

  for (seed in 1:run) {
    print(seed)
    data.sim <- main.fun(I, J, n, sam.size, tau3, tau2, sigma2, alpha.0, alpha.1, alpha.2, psi.0, psi.1, psi.2, psi.3, seed)

    #Generate testing dataset
    data.test <- main.fun(I, J, n, sam.size, tau3, tau2, sigma2, alpha.0, alpha.1, alpha.2, psi.0, psi.1, psi.2, psi.3, seed+1)

    #Extract summary information on number of alters in the test data
    Min.Num.Alters <- c(Min.Num.Alters, min(data.test$num.alters))
    Max.Num.Alters <- c(Max.Num.Alters, max(data.test$num.alters))
    Mean.Num.Alters <- c(Mean.Num.Alters, mean(data.test$num.alters))
    Median.Num.Alters <- c(Median.Num.Alters, median(data.test$num.alters))
    
    A1opt <- data.test$AOpt.ijk.1
    A2opt <- data.test$AOpt.ijk.2
    
    Aopt.dwols.both <- dwols.both(data.sim, data.test, I, J, n)
    A1opt.dwols.both <- Aopt.dwols.both$Aopt.dwols.1.test
    A2opt.dwols.both <- Aopt.dwols.both$Aopt.dwols.2.test
   
    Yopt.dwols.both <- Aopt.dwols.both$Yopt.test
    
    Yobs <- data.sim$Y.ijk
    
    Yopt <- data.test$Y.ijk + (A1opt - data.test$A.ijk.1)*(psi.0 + psi.1*data.test$X.ijk.1 + psi.2*data.test$total.t.1 + psi.3*data.test$V.i.1*data.test$total.t.1) + (A2opt - data.test$A.ijk.2)*(psi.0 + psi.1*data.test$X.ijk.2 + psi.2*data.test$total.t.2 + psi.3*data.test$V.i.2*data.test$total.t.2)

    Aopt.dwols.both.resc <- dwols.both.resc(data.sim, data.test, I, J, n)
    A1opt.dwols.both.resc <- Aopt.dwols.both.resc$Aopt.dwols.1.test
    A2opt.dwols.both.resc <- Aopt.dwols.both.resc$Aopt.dwols.2.test
   
    Yopt.dwols.both.resc <- Aopt.dwols.both.resc$Yopt.test
    
    Yobs <- data.sim$Y.ijk
    
    Yopt <- data.test$Y.ijk + (A1opt - data.test$A.ijk.1)*(psi.0 + psi.1*data.test$X.ijk.1 + psi.2*data.test$total.t.1 + psi.3*data.test$V.i.1*data.test$total.t.1) + (A2opt - data.test$A.ijk.2)*(psi.0 + psi.1*data.test$X.ijk.2 + psi.2*data.test$total.t.2 + psi.3*data.test$V.i.2*data.test$total.t.2)

    Aopt.dwols.wthn <- dwols.wthn(data.sim, data.test, I, J, n)
    A1opt.dwols.wthn <- Aopt.dwols.wthn$Aopt.dwols.1.test
    A2opt.dwols.wthn <- Aopt.dwols.wthn$Aopt.dwols.2.test
   
    Yopt.dwols.wthn <- Aopt.dwols.wthn$Yopt.test

    Aopt.dwols.standard <- dwols.standard(data.sim, data.test, I, J, n)
    A1opt.dwols.standard<- Aopt.dwols.standard$Aopt.dwols.1.test
    A2opt.dwols.standard <- Aopt.dwols.standard$Aopt.dwols.2.test
   
    Yopt.dwols.standard <- Aopt.dwols.standard$Yopt.test
    
    Res.both[seed,] <- c(c(sum(data.sim$AOpt.ijk.1 == data.sim$A.ijk.1), sum(data.sim$AOpt.ijk.2 == data.sim$A.ijk.2), sum(A1opt.dwols.both == A1opt), sum(A2opt.dwols.both == A2opt))/(sam.size*I), mean(Yopt), mean(Yopt.dwols.both), Aopt.dwols.both$psi.0.stage1[1], Aopt.dwols.both$psi.1.stage1[1], Aopt.dwols.both$psi.2.stage1[1], Aopt.dwols.both$psi.3.stage1[1], Aopt.dwols.both$psi.0.stage2[1], Aopt.dwols.both$psi.1.stage2[1], Aopt.dwols.both$psi.2.stage2[1], Aopt.dwols.both$psi.3.stage2[1], Aopt.dwols.both$xsi.0.stage1[1], Aopt.dwols.both$xsi.1.stage1[1], Aopt.dwols.both$xsi.0.stage2[1], Aopt.dwols.both$xsi.1.stage2[1])
    Res.both.resc[seed,] <- c(c(sum(data.sim$AOpt.ijk.1 == data.sim$A.ijk.1), sum(data.sim$AOpt.ijk.2 == data.sim$A.ijk.2), sum(A1opt.dwols.both.resc == A1opt), sum(A2opt.dwols.both.resc == A2opt))/(sam.size*I), mean(Yopt), mean(Yopt.dwols.both.resc), Aopt.dwols.both.resc$psi.0.stage1[1], Aopt.dwols.both.resc$psi.1.stage1[1], Aopt.dwols.both.resc$psi.2.stage1[1], Aopt.dwols.both.resc$psi.3.stage1[1], Aopt.dwols.both.resc$psi.0.stage2[1], Aopt.dwols.both.resc$psi.1.stage2[1], Aopt.dwols.both.resc$psi.2.stage2[1], Aopt.dwols.both.resc$psi.3.stage2[1], Aopt.dwols.both.resc$xsi.0.stage1[1], Aopt.dwols.both.resc$xsi.1.stage1[1], Aopt.dwols.both.resc$xsi.0.stage2[1], Aopt.dwols.both.resc$xsi.1.stage2[1])
    Res.wthn[seed,] <- c(c(sum(data.sim$AOpt.ijk.1 == data.sim$A.ijk.1), sum(data.sim$AOpt.ijk.2 == data.sim$A.ijk.2), sum(A1opt.dwols.wthn == A1opt), sum(A2opt.dwols.wthn == A2opt))/(sam.size*I), mean(Yopt), mean(Yopt.dwols.wthn), Aopt.dwols.wthn$psi.0.stage1[1], Aopt.dwols.wthn$psi.1.stage1[1], Aopt.dwols.wthn$psi.2.stage1[1], Aopt.dwols.wthn$psi.3.stage1[1], Aopt.dwols.wthn$psi.0.stage2[1], Aopt.dwols.wthn$psi.1.stage2[1], Aopt.dwols.wthn$psi.2.stage2[1], Aopt.dwols.wthn$psi.3.stage2[1], Aopt.dwols.wthn$xsi.0.stage1[1], Aopt.dwols.wthn$xsi.1.stage1[1], Aopt.dwols.wthn$xsi.0.stage2[1], Aopt.dwols.wthn$xsi.1.stage2[1])
    Res.standard[seed,] <- c(c(sum(data.sim$AOpt.ijk.1 == data.sim$A.ijk.1), sum(data.sim$AOpt.ijk.2 == data.sim$A.ijk.2), sum(A1opt.dwols.standard == A1opt), sum(A2opt.dwols.standard == A2opt))/(sam.size*I), mean(Yopt), mean(Yopt.dwols.standard), Aopt.dwols.standard$psi.0.stage1[1], Aopt.dwols.standard$psi.1.stage1[1], Aopt.dwols.standard$psi.2.stage1[1], Aopt.dwols.standard$psi.3.stage1[1], Aopt.dwols.standard$psi.0.stage2[1], Aopt.dwols.standard$psi.1.stage2[1], Aopt.dwols.standard$psi.2.stage2[1], Aopt.dwols.standard$psi.3.stage2[1], Aopt.dwols.standard$xsi.0.stage1[1], Aopt.dwols.standard$xsi.1.stage1[1], Aopt.dwols.standard$xsi.0.stage2[1], Aopt.dwols.standard$xsi.1.stage2[1])

  }
  return(list("Res.both"=Res.both, "Res.both.resc"=Res.both.resc, "Res.wthn"=Res.wthn, "Res.standard"=Res.standard, "Min.Num.Alters" = Min.Num.Alters, "Max.Num.Alters" = Max.Num.Alters, "Mean.Num.Alters" = Mean.Num.Alters, "Median.Num.Alters" = Median.Num.Alters))
}





#Only need to iterate through each stage once
p.list <- c(1)

#Storage
A.Opts.1 <- matrix(nrow = length(p.list), ncol = 12)
A.Opts.2 <- matrix(nrow = length(p.list), ncol = 12)
A.Opts.3 <- matrix(nrow = length(p.list), ncol = 12)
A.Opts.4 <- matrix(nrow = length(p.list), ncol = 12)

Y.Opts.median.1 <- matrix(nrow = length(p.list), ncol = 3)
Y.Opts.median.2 <- matrix(nrow = length(p.list), ncol = 3)
Y.Opts.median.3 <- matrix(nrow = length(p.list), ncol = 3)
Y.Opts.median.4 <- matrix(nrow = length(p.list), ncol = 3)

#Counter to mark progress
count <- 0

run <- 1000
I <- 10
J <- 2
n <- 50
sam.size <- 50

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
#Only need to iterate through each stage once
for (p in p.list) {
  count <- count + 1
  #Run simulations
  Results <- Sims(I, J, n, sam.size, tau3, tau2, sigma2, alpha.0, alpha.1, alpha.2, psi.0, psi.1, psi.2, psi.3, run, model)
  Results.both <- Results$Res.both
  Results.both.resc <- Results$Res.both.resc
  Results.wthn <- Results$Res.wthn
  Results.standard <- Results$Res.standard

  Min.Num.Alters <- Results$Min.Num.Alters
  Max.Num.Alters <- Results$Max.Num.Alters
  Mean.Num.Alters <- Results$Mean.Num.Alters
  Median.Num.Alters <- Results$Median.Num.Alters

  cat('count ', count, ' done\n')
}

outputs <- list(
  TwoStage_results_both_trt_fr_mis <- Results.both,
  TwoStage_results_both__resc_trt_fr_mis <- Results.both.resc,
  TwoStage_results_wthn_trt_fr_mis <- Results.wthn,
  TwoStage_results_standard_trt_fr_mis <- Results.standard,
  TwoStage_min_num_trt_fr_mis <- Min.Num.Alters,
  TwoStage_max_num_trt_fr_mis <- Max.Num.Alters,
  TwoStage_mean_num_trt_fr_mis <- Mean.Num.Alters,
  TwoStage_median_num_trt_fr_mis <- Median.Num.Alters
)

saveRDS(outputs, "results/sims_two_stage_trt_fr_mis.rds")

dwols.trt.df <- data.frame(AOpts = 100*c(c(Results.wthn$DWOLS.WTHN.STAGE1),c(Results.wthn$DWOLS.WTHN.STAGE2),c(Results.standard$DWOLS.STANDARD.STAGE1),c(Results.standard$DWOLS.STANDARD.STAGE2),c(Results.both$DWOLS.BOTH.STAGE1),c(Results.both$DWOLS.BOTH.STAGE2),c(Results.both.resc$DWOLS.BOTH.RESC.STAGE1),c(Results.both.resc$DWOLS.BOTH.RESC.STAGE2)), Method = rep(c("ME1","ME2","ME3", "ME4"),each=run*2), Stage = rep(rep(c("Stage 1","Stage 2"),each=run),4) )

ggplot(dwols.trt.df, aes(x=Stage, y=AOpts, fill=Method)) + geom_boxplot() + scale_fill_grey(start=0.8, end=0.4) + xlab("Treatment Stage") + ylab("Optimal Treatment Rate (%)")
#Plot outcome results

dwols.out.df <- data.frame(YOpts = c(c(Results.both$Yopt), c(Results.wthn$Yopt.DWOLS.WTHN), c(Results.standard$Yopt.DWOLS.STANDARD), c(Results.both$Yopt.DWOLS.BOTH),c(Results.both.resc$Yopt.DWOLS.BOTH.RESC)), Outcome = rep(c("True", "ME1","ME2", "ME3","ME4"),each=1000))

dwols.out.df$Outcome <- factor(dwols.out.df$Outcome, levels=c("ME1","ME2", "ME3", "ME4","True"))
ggplot(dwols.out.df, aes(y=YOpts, fill=Outcome)) + 
  geom_boxplot() + scale_fill_grey(start=0.8, end=0.4) + ylab("Optimal Outcome")


dwols.psi.stage1.df <- data.frame(Psi = c(c(Results.both.resc$PSI.0.STAGE1.BOTH.RESC),c(Results.both.resc$PSI.1.STAGE1.BOTH.RESC),c(Results.both.resc$PSI.2.STAGE1.BOTH.RESC),c(Results.both.resc$PSI.3.STAGE1.BOTH.RESC)), Blip = rep(c("Psi0", "Psi1", "Psi2", "Psi3"),each=1000))
dwols.psi.stage2.df <- data.frame(Psi = c(c(Results.both.resc$PSI.0.STAGE2.BOTH.RESC),c(Results.both.resc$PSI.1.STAGE2.BOTH.RESC),c(Results.both.resc$PSI.2.STAGE2.BOTH.RESC),c(Results.both.resc$PSI.3.STAGE2.BOTH.RESC)), Blip = rep(c("Psi0", "Psi1", "Psi2", "Psi3"),each=1000))

dwols.psi.stage1.df$Blip <- factor(dwols.psi.stage1.df$Blip, levels=c("Psi0", "Psi1", "Psi2", "Psi3"))
dwols.psi.stage2.df$Blip <- factor(dwols.psi.stage2.df$Blip, levels=c("Psi0", "Psi1", "Psi2", "Psi3"))

ggplot(dwols.psi.stage1.df, aes(y=Psi, fill=Blip)) + 
  geom_boxplot() + scale_fill_grey(start=0.8, end=0.4) + ylab("Blip Parameters: Stage 1")

ggplot(dwols.psi.stage2.df, aes(y=Psi, fill=Blip)) + 
  geom_boxplot() + scale_fill_grey(start=0.8, end=0.4) + ylab("Blip Parameters: Stage 2")


standard.psi.stage1.df <- data.frame(Psi = c(c(Results.standard$PSI.0.STAGE1.STANDARD),c(Results.standard$PSI.1.STAGE1.STANDARD),c(Results.standard$PSI.2.STAGE1.STANDARD),c(Results.standard$PSI.3.STAGE1.STANDARD)), Blip = rep(c("Psi0", "Psi1", "Psi2", "Psi3"),each=1000))
standard.psi.stage2.df <- data.frame(Psi = c(c(Results.standard$PSI.0.STAGE2.STANDARD),c(Results.standard$PSI.1.STAGE2.STANDARD),c(Results.standard$PSI.2.STAGE2.STANDARD),c(Results.standard$PSI.3.STAGE2.STANDARD)), Blip = rep(c("Psi0", "Psi1", "Psi2", "Psi3"),each=1000))

standard.psi.stage1.df$Blip <- factor(standard.psi.stage1.df$Blip, levels=c("Psi0", "Psi1", "Psi2", "Psi3"))
standard.psi.stage2.df$Blip <- factor(standard.psi.stage2.df$Blip, levels=c("Psi0", "Psi1", "Psi2", "Psi3"))

ggplot(standard.psi.stage1.df, aes(y=Psi, fill=Blip)) + 
  geom_boxplot() + scale_fill_grey(start=0.8, end=0.4) + ylab("Blip Parameters: Stage 1")

ggplot(standard.psi.stage2.df, aes(y=Psi, fill=Blip)) + 
  geom_boxplot() + scale_fill_grey(start=0.8, end=0.4) + ylab("Blip Parameters: Stage 2")






#***************************************************************************************************
