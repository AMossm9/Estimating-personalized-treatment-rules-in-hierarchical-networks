#**********DWOLS FOR WITHIN- AND BETWEEN-GROUP INTERFERENCE: ONE-STAGE SETTING, TREATMENT AND TREATMENT-FREE MODELS MISSPECIFIED**********

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
  trt.df <- matrix(, nrow=N, ncol=8)
  index <- 1

  #Iterate through each hierarchy (i.e. each Level 3 group)
  for(i in 1:I){

    #Generate random intercept for Level 3
    beta.0i <- sqrt(tau3)*rnorm(1, 0, 1)

    #Generate cluster-level covariate for Level 3
    V.i <- rnorm(1, 1, 1)

    #Iterate through each Level 2 group
    for(j in 1:J){

      #Generate random intercept for Level 2
      beta.0j <- sqrt(tau2)*rnorm(1, 0, 1)

      #Iterate through each unit in Level 2
      for(k in 1:n){

        #Generate an individual level covariate
        X.ijk <- rnorm(1, 2, 1)

        #Generate error for treatment model
        eps.ijk <- sqrt(sigma2)*rnorm(1, 0, 1)

        #Treatment model with random intercepts: this is the TRUE UNDERLYING TREATMENT MODEL that we will incorrectly specify
        #in our estimation
        alpha.3 <- 0.5*alpha.2
        logit.p <- alpha.1*X.ijk + alpha.2*V.i + alpha.3*(sin(X.ijk)) + alpha.0 + beta.0j + beta.0i

        #Treatment model without random intercepts (hence "fixed")
        logit.p.fixed <- alpha.1*X.ijk + alpha.2*V.i + alpha.3*(sin(X.ijk)) + alpha.0

        #Probability of treatment (trt) is based on true underlying treatment model
        prob.trt <- exp(logit.p)/(exp(logit.p)+1)
        A.ijk <- rbinom(1, 1, prob.trt)

        #Indicators for whether an individual receives treatment or not
        #Not sure if we need this later on
        trt.indicator <- ifelse(A.ijk == 1, 1, -1)

        #Store data
        trt.df[index, 1:8] <- c(beta.0i, V.i, beta.0j, X.ijk, trt.indicator, logit.p.fixed, logit.p, A.ijk)
        
        index <- index + 1
      }
    }
  }

  trt.data <- cbind(L3, L2, trt.df)
  trt.data <- data.frame(trt.data)
  #Generate the unique N different patient ID's
  trt.data$ID <- seq(1, N, 1)

  names(trt.data) <- c("L3", "L2", "beta.0i", "V.i", "beta.0j", "X.ijk", "trt.indicator", "logit.p.fixed", "logit.p", "A.ijk", "ID")

  #MIS-SPECIFIED TREATMENT MODEL IGNORES NON-LINEAR TERM IN X.ijk
  #Fixed-effects only model
  mod.1 <- glm(A.ijk~X.ijk + V.i, family = 'binomial', data = trt.data)
  
  #Fixed- and random-effects model 
  mod.2 <- glmer(A.ijk~X.ijk + V.i + (1 | L3/L2), family = 'binomial', data = trt.data)

  #summary(mod.2)
 
  #Store the estimated standard deviations of the random effects
  #Assuming independence of the random effects, we can sum the standard deviations for the standard deviation of 
  #the normal distribution of the random effects
  rand.sd <- sum(as.data.frame(VarCorr(mod.2))[1:2,5])

  #IMPORTANT: now, logit.p.fixed and logit.p are different from what we defined earlier
  #logit.p.fixed is the ESTIMATED logit of probability of treatment, integrating out the random effects of mod.2
  #Refer to stackexchange for why we set re.form=NA
  #Integrating out the random effects seems to be consistent with how the joint propensity score with random effects
  #is defined in the Perez-Heydrich et al. (2014) paper
  #logit.p is the ESTIMATED logit if probability of treatment using mod.1
  trt.data$logit.p.fixed <- predict(mod.2, trt.data, type='response', re.form=NULL)
  
  trt.data$logit.p <- predict(mod.1, trt.data, type='response')
  
  #Return a dataframe with N x 11 columns 
  return(list(trt.data=trt.data, rand.sd=rand.sd))
}

#**********EGO AND ALTER SELECTION**********

#Takes seed, generated treatment data, number of Level 3 groups, number of Level 2 groups,
#number of units in each Level 2 group, and number of egos (sam.size) we want to select

egos.alters <- function(seed, trt.data, I, J, n, sam.size){

  set.seed(seed)

  #List to store data for the egos selected
  #List to store data for the remaining alters
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
    
    #Randomly assign alters to egos so that EACH EGO has AT LEAST ONE ALTER
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
    df.alters$alter <- rep("alter", nrow(df.alters))
    
    mydata.egos[[i]] <- df.egos.rep
    mydata.alters[[i]] <- df.alters
    
  }

  data.egos <- do.call(rbind, mydata.egos)
  #Distinguish which L2 group the egos belong to
  names(data.egos)[2] <- "ego.L2"
  data.alters <- do.call(rbind, mydata.alters)
  #Distinguish which L2 group the alters belong to
  names(data.alters)[2] <- "alter.L2"
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
  #Rename columns 11 and 23 so that we're distinguishing which ID's correspond to
  #egos and which ID's correspond to alters
  names(trt.data)[11] <- "ID.Ego"
  names(trt.data)[23] <- "ID.Alter"

  #alters.joint() was initially a product of alters' propensities, later becomes total number of treated
  #alters (accounting for both within- and between-group status)
  #For generating treatment decision rules and outcomes, we want to keep 
  #track of total number of within-group and between-group neighbours an ego has.
  alters.joint <- c()
  total.w <- c()
  total.b <- c()

  #Iterate through egos, first through Level 3 and then through Level 2 in Level 3
  for(i in 1:I){

    df.L3 <- trt.data[trt.data$L3 == i,]
    ego.w <- c()
    ego.b <- c()

    #This will be used for alters.joint().
    alters.prop.prod <- c()

    #We need to search through each unique ego ID in our dataframe (since
    #there are multiple rows for an ego, each linked with a different alter's information).
    unique.egos <- unique(df.L3$ID.Ego)
    for(u in unique.egos){

      #Extract dataframe for that ego
      df.ego <- df.L3[df.L3$ID.Ego==u,]

      #Find its within-group neighbours
      df.ego.w <- df.ego[df.ego$alter.group == "W",]

      #Sum of treatment status of within-group neighbours
      sum.trt.w <- sum(df.ego.w[,22])

      #Find its between-group neighbours
      df.ego.b <- df.ego[df.ego$alter.group == "B",]

      #Sum of treatment status of between-group neighbours
      sum.trt.b <- sum(df.ego.b[,22])

      ego.w <- c(ego.w, rep(sum.trt.w, nrow(df.ego)))

      ego.b <- c(ego.b, rep(sum.trt.b, nrow(df.ego)))

      #Take the product of the ego's alters' propensities
      alters.prop.prod <- c(alters.prop.prod, rep(prod(df.ego[,21]), nrow(df.ego)))
    }

  alters.joint <- c(alters.joint, alters.prop.prod)
  total.w <- c(total.w, ego.w)
  total.b <- c(total.b, ego.b)
  }
  
  trt.data$alters.joint <- alters.joint
  trt.data$total.w <- total.w
  trt.data$total.b <- total.b
  #This is where alters.joint turns into the total number of treated neighbours,
  #denoted with "total.t"
  trt.data$alters.joint <- total.w + total.b
  colnames(trt.data)[which(names(trt.data) == "alters.joint")] <- "total.t"

  #Initialize empty vectors for storing outcome information 
  Y.ijk <- c()
  Y.err.ijk <- c()
  AOpt.ijk <- c()
  tf.ego <- c()
  blip.ego <- c()

  #Iterate through Level 3 groups
  for(i in 1:I){

    #Extract ego-specific information needed
    df.L3 <- trt.data[trt.data$L3 == i,c(1:12,23:28)]
    Y <- c()
    AOpt <- c()
    Y.err <- c()
    blips <- c()
    tf <- c()
    #Iterate through unique ego ID's
    unique.egos <- unique(df.L3$ID.Ego)
    for(u in unique.egos){
      #Extract a given ego's information (all of their rows, this is fine
      #since their cluster- and patient-specific information are the same across rows).
      df.ego <- df.L3[df.L3$ID.Ego==u,]
      #print(df.ego)

      #Treatment-free model: comprised of individual and cluster-level covariates
      #Interference terms account for both within- and between-groups of neighbours (i.e. total.t)
      #Interference term is interacting with cluster-level covariate V.i
      trt.free <- df.ego$X.ijk + (df.ego$X.ijk)^2 + df.ego$V.i + df.ego$total.t + df.ego$V.i*df.ego$total.t

      #Optimal treatment rule based on blip component
      A.Opt <- as.numeric((psi.0 + psi.1*df.ego$X.ijk + psi.2*df.ego$total.t + psi.3*df.ego$V.i*df.ego$total.t)>0)

      #Calculate regrets and blips
      mu <- (A.Opt-df.ego$A.ijk)*(psi.0 + psi.1*df.ego$X.ijk + psi.2*df.ego$total.t + psi.3*df.ego$V.i*df.ego$total.t)

      blip <- (df.ego$A.ijk)*(psi.0 + psi.1*df.ego$X.ijk + psi.2*df.ego$total.t + psi.3*df.ego$V.i*df.ego$total.t)

      #Overall outcome model
      mu.ijk <- trt.free + blip

      #Error component for outcome model
      y.err <- rep(rnorm(1, 0, 1), length(mu.ijk))
      #Y.out <- rnorm(length(mu.ijk), mu.ijk, 1)
      Y.out <- rep(rnorm(1, mu.ijk[1], 1), length(mu.ijk))

      Y <- c(Y, Y.out)
      AOpt <- c(AOpt, A.Opt)
      Y.err <- c(Y.err, y.err)
      blips <- c(blips, blip)
      tf <- c(tf, trt.free)
    }

    Y.ijk <- c(Y.ijk, Y)
    AOpt.ijk <- c(AOpt.ijk, AOpt)
    Y.err.ijk <- c(Y.err.ijk, Y.err)
    blip.ego <- c(blip.ego, blips)
    tf.ego <- c(tf.ego, tf)
  }

  trt.data$Y.ijk <- Y.ijk
  trt.data$AOpt.ijk <- AOpt.ijk
  trt.data$Y.err.ijk <- Y.err.ijk
  trt.data$blip.ego <- blip.ego
  trt.data$tf.ego <- tf.ego
  #print(head(trt.data))
  return(trt.data)
}

#**********WEIGHT CONSTRUCTION**********

#The group.propensity function will compute different weights 
#used in our different methods

group.propensity <- function(trt.data, I, J, n, tau3, tau2, rand.sd){

  #w.L3: network absolute value weights that account for both within- and between-group neighbours
  w.L3 <- c()
  #w.within: network absolute value weights that account for only within-group neighbours
  w.within <- c()
  #w.standard: calculates the ego's absolute value weight using fixed-effects only model,
  #does not account for neighbours 
  w.standard <- c()
  #w.within: network absolute value weights that account for only within-group neighbours
  w.network <- c()

  #Subset data based on Level 3
  for(i in 1:I){
    df.L3 <- trt.data[trt.data$L3 == i,]
    #Use ".L2" as we start storing weights
    #because we are now working with 
    #egos in Level 2 groups within a given
    #Level 3 group
    w.L2 <- c()
    w.within.L2 <- c()
    w.standard.L2 <- c()
    w.network.L2 <- c()

    unique.egos <- unique(df.L3$ID.Ego)
    
    #Iterate through each ego within Level 3 (possibly
    #from different Level 2 groups)
    for(u in unique.egos){

      #Get dataframe corresponding to ego
      df.ego <- df.L3[df.L3$ID.Ego==u,]

      #Get the ego's probability of treatment from mixed model (it's called logit but is actually the expit)
      #and indicator for whether or not they received treatment "1"
      ego.logit <- df.ego[,8][1]
      ego.indicator <- df.ego[,7][1]
      ego.A <- ifelse(ego.indicator==1, 1, 0)

      #Get vector of alters' (within- and between-group) probability of treatment from mixed model
      alter.logit <- df.ego[,20]
      alter.indicator <- df.ego[,19]
      alter.A <- ifelse(alter.indicator==1, 1, 0)

      #Calculate network absolute value Weights for group propensity score
      w.L2 <- c(w.L2, rep(abs(ego.A-ego.logit)*prod(abs(alter.A-alter.logit)), nrow(df.ego)))

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
        ego.w.logit <- df.ego[,8][1]
        ego.w.indicator <- df.ego[,7][1]
        ego.w.A <- ifelse(ego.w.indicator==1, 1, 0)
      
        w.within.L2 <- c(w.within.L2, rep(abs(ego.w.A-ego.w.logit), nrow(df.ego)))
        w.network.L2 <- c(w.network.L2, rep(abs(ego.w.A-ego.w.logit), nrow(df.ego)))

      }

      else{
        ego.w.logit <- df.ego[,8][1]
        ego.w.indicator <- df.ego[,7][1]
        ego.w.A <- ifelse(ego.w.indicator==1, 1, 0)
        alter.w.logit <- df.ego.w[,20]
        alter.w.indicator <- df.ego.w[,19]
        alter.w.A <- ifelse(alter.w.indicator==1, 1, 0)
        
        w.within.L2 <- c(w.within.L2, rep(abs(ego.w.A-ego.w.logit)*prod(abs(alter.w.A-alter.w.logit)), nrow(df.ego)))
        w.network.L2 <- c(w.network.L2, rep(abs(ego.w.A-ego.w.logit)*prod(abs(alter.w.A-alter.w.logit)), nrow(df.ego)))
      }

      #Weights for standard individual absolute value weight
      #from mixed-effects model
      ego.s.logit <- df.ego[,8][1]
      ego.s.indicator <- df.ego[,7][1]
      ego.s.A <- ifelse(ego.s.indicator==1, 1, 0)
      ipw.w <- abs(ego.s.A - ego.s.logit)
      w.standard.L2 <- c(w.standard.L2, rep(ipw.w, nrow(df.ego)))
    }

    w.L3 <- c(w.L3, w.L2)
    w.within <- c(w.within, w.within.L2)
    w.standard <- c(w.standard, w.standard.L2)
    w.network <- c(w.network, w.network.L2)
  }

  trt.data$w.L3 <- w.L3
  trt.data$w.within <- w.within
  trt.data$w.standard <- w.standard
  trt.data$w.network <- w.network
  return(trt.data)
}

#**********GENERATE DATASET FOR EGOS**********

#Initializes hierarchical network structure, selects egos and alters,
#generates outcomes, and calculates different weights
main.fun <- function(I, J, n, sam.size, tau3, tau2, sigma2, alpha.0, alpha.1, alpha.2, psi.0, psi.1, psi.2, psi.3, seed){
  data.gen.list <- data.gen(seed, I, J, n, tau3, tau2, sigma2, alpha.0, alpha.1, alpha.2)
  trt.data <- data.gen.list[[1]]
  rand.sd <- data.gen.list[[2]]
  trt.data <- egos.alters(seed, trt.data, I, J, n, sam.size)
  trt.data <- egos.outcomes(trt.data, I, J, n, psi.0, psi.1, psi.2, psi.3)
  trt.data <- group.propensity(trt.data, I, J, n, tau3, tau2, rand.sd)
  
  ego.df <- trt.data[, c(1:11, 26:37)]

  num.reps <- as.numeric(table(ego.df$ID.Ego))

  ego.df <- ego.df[order(ego.df$ID.Ego),]

  ego.df$num.alters <- rep(num.reps, num.reps) 

  ego.df <- ego.df[!duplicated(ego.df$ID.Ego), ]

  ego.df <- ego.df[order(ego.df$num.alters),]

  w.within.mean <- c()

  w.network.mean <- c()

  w.L3.mean <- c()

  for(t in unique(ego.df$num.alters)){

      df.t <- ego.df[which(ego.df$num.alters==t),]

      w.within.mean <- c(w.within.mean, rep(mean(df.t$w.within), nrow(df.t)))
      
      w.network.mean <- c(w.network.mean, rep(mean(df.t$w.network), nrow(df.t)))
      
      w.L3.mean <- c(w.L3.mean, rep(mean(df.t$w.L3), nrow(df.t)))
      
  }
  
  ego.df$w.within.mean <- w.within.mean
  ego.df$w.network.mean <- w.network.mean
  ego.df$w.L3.mean <- w.L3.mean

  ego.df$w.within.resc <- ego.df$w.within/w.within.mean

  ego.df$w.network.resc <- ego.df$w.network/w.network.mean

  ego.df$w.L3.resc <- ego.df$w.L3/w.L3.mean
  
  return(ego.df)
}

#Test out a few different datasets
ego.df1 <- main.fun(10, 4, 20, 20, 1, 1, 1, -1, 1, 0.5, 1, 1, -1, -1.5, 1000)
sum(abs(ego.df1$AOpt.ijk-ego.df1$A.ijk))/nrow(ego.df1)

ego.df2 <- main.fun(10, 5, 20, 20, 1, 1, 1, -1, 1, 0.5, 1, 1, -1, -1.5, 1000)
sum(abs(ego.df2$AOpt.ijk-ego.df2$A.ijk))/nrow(ego.df2)

ego.df3 <- main.fun(10, 6, 20, 20, 1, 1, 1, -1, 1, 0.5, 1, 1, -1, -1.5, 1000)
sum(abs(ego.df3$AOpt.ijk-ego.df3$A.ijk))/nrow(ego.df3)

#**********DWOLS IMPLEMENTATION WHERE TREATMENT-FREE MODEL IS MIS-SPECIFIED VIA OMISSION OF NON-LINEAR TERM IN X.ijk**********

#METHOD #1: ONLY ACCOUNTS FOR WITHIN-GROUP INTERFERENCE, USES WITHIN-GROUP WEIGHTS

dwols <- function(data.sim, data.test, I, J, n) {
  N <- I*J*n

  mod <- lm(Y.ijk ~ X.ijk + V.i + total.w + V.i:total.w + A.ijk + A.ijk:X.ijk + A.ijk:total.w + A.ijk:V.i:total.w, weights = w.within, data = data.sim)
  
  #Obtain extra blip parameter estimates from outcome model
  #(only coefficients for treatment terms)
  psi.0.hat <- coef(mod)[5]
  psi.1.hat <- coef(mod)[7]
  psi.2.hat <- coef(mod)[8]
  psi.3.hat <- coef(mod)[9]

  #Obtain spillover parameter estimates in the treatment-free model
  xsi.0.hat <- coef(mod)[4]
  xsi.1.hat <- coef(mod)[6]

  #Compute optimal treatment in testing set using blip parameter estimates
  Aopt.dwols <- as.numeric((psi.0.hat + psi.1.hat*data.test$X.ijk + psi.2.hat*data.test$total.w + psi.3.hat*data.test$V.i*data.test$total.w) > 0)
  
  #Optimal outcome
  Y.opt <- data.test$Y.ijk + (Aopt.dwols - data.test$A.ijk)*(psi.0.hat + psi.1.hat*data.test$X.ijk + psi.2.hat*data.test$total.w + psi.3.hat*data.test$V.i*data.test$total.w)
  
  #Regret
  regret <- (Aopt.dwols - data.test$A.ijk)*(psi.0.hat + psi.1.hat*data.test$X.ijk + psi.2.hat*data.test$total.w + psi.3.hat*data.test$V.i*data.test$total.w)

  #E[Y|A = Aopt.dwols]
  #Optimal treatment based on what we observe in given dataset
  Aopt.dwols.obs <- as.numeric((psi.0.hat + psi.1.hat*data.sim$X.ijk + psi.2.hat*data.sim$total.w + psi.3.hat*data.sim$V.i*data.sim$total.w) > 0)
  
  #Mean of outcomes where treatment we observe is the same as the dwols optimal
  #treatment for our given dataset
  Y.agree <- mean(data.sim$Y.ijk[data.sim$A.ijk == Aopt.dwols.obs])
  
  return(data.frame(mean(Aopt.dwols == data.test$AOpt.ijk), mean(Y.opt), mean(Y.agree), psi.0.hat, psi.1.hat, psi.2.hat, psi.3.hat, xsi.0.hat, xsi.1.hat, mean(regret)))
}

#METHOD #2: ACCOUNTS FOR WITHIN- AND BETWEEN-GROUP INTERFERENCE, USES STANDARD INDIVIDUAL WEIGHTS

dwols_mod_1 <- function(data.sim, data.test, I, J, n) {
  N <- I*J*n

  mod <- lm(Y.ijk ~ X.ijk + V.i + total.t + V.i:total.t + A.ijk + A.ijk:X.ijk + A.ijk:total.t + A.ijk:V.i:total.t, weights = w.standard, data = data.sim)
  
  #Obtain extra blip parameter estimates from outcome model
  #(only coefficients for treatment terms)
  psi.0.hat <- coef(mod)[5]
  psi.1.hat <- coef(mod)[7]
  psi.2.hat <- coef(mod)[8]
  psi.3.hat <- coef(mod)[9]

  #Obtain spillover parameter estimates in the treatment-free model
  xsi.0.hat <- coef(mod)[4]
  xsi.1.hat <- coef(mod)[6]

  #Compute optimal treatment in testing set using blip parameter estimates
  Aopt.dwols <- as.numeric((psi.0.hat + psi.1.hat*data.test$X.ijk + psi.2.hat*data.test$total.t + psi.3.hat*data.test$V.i*data.test$total.t) > 0)
  
  #Optimal outcome
  Y.opt <- data.test$Y.ijk + (Aopt.dwols - data.test$A.ijk)*(psi.0.hat + psi.1.hat*data.test$X.ijk + psi.2.hat*data.test$total.t + psi.3.hat*data.test$V.i*data.test$total.t)
  
  #Regret
  regret <- (Aopt.dwols - data.test$A.ijk)*(psi.0.hat + psi.1.hat*data.test$X.ijk + psi.2.hat*data.test$total.t + psi.3.hat*data.test$V.i*data.test$total.t)

  #E[Y|A = Aopt.dwols]
  #Optimal treatment based on what we observe in given dataset
  Aopt.dwols.obs <- as.numeric((psi.0.hat + psi.1.hat*data.sim$X.ijk + psi.2.hat*data.sim$total.t + psi.3.hat*data.sim$V.i*data.sim$total.t) > 0)
  
  #Mean of outcomes where treatment we observe is the same as the dwols optimal
  #treatment for our given dataset
  Y.agree <- mean(data.sim$Y.ijk[data.sim$A.ijk == Aopt.dwols.obs])
  
  return(data.frame(mean(Aopt.dwols == data.test$AOpt.ijk), mean(Y.opt), mean(Y.agree), psi.0.hat, psi.1.hat, psi.2.hat, psi.3.hat, xsi.0.hat, xsi.1.hat, mean(regret)))
}

#METHOD #3: ACCOUNTS FOR WITHIN- AND BETWEEN-GROUP INTERFERENCE, USES NON-RESCALED NETWORK WEIGHTS FROM ALL NEIGHBOURS

dwols_mod_2 <- function(data.sim, data.test, I, J, n) {
  N <- I*J*n

  mod <- lm(Y.ijk ~ X.ijk + V.i + total.t + V.i:total.t + A.ijk + A.ijk:X.ijk + A.ijk:total.t + A.ijk:V.i:total.t, weights = w.L3, data = data.sim)
  
  #Obtain extra blip parameter estimates from outcome model
  #(only coefficients for treatment terms)
  psi.0.hat <- coef(mod)[5]
  psi.1.hat <- coef(mod)[7]
  psi.2.hat <- coef(mod)[8]
  psi.3.hat <- coef(mod)[9]

  #Obtain spillover parameter estimates in the treatment-free model
  xsi.0.hat <- coef(mod)[4]
  xsi.1.hat <- coef(mod)[6]
  
  #Compute optimal treatment in testing set using blip parameter estimates
  Aopt.dwols <- as.numeric((psi.0.hat + psi.1.hat*data.test$X.ijk + psi.2.hat*data.test$total.t + psi.3.hat*data.test$V.i*data.test$total.t) > 0)
  
  #Optimal outcome
  Y.opt <- data.test$Y.ijk + (Aopt.dwols - data.test$A.ijk)*(psi.0.hat + psi.1.hat*data.test$X.ijk + psi.2.hat*data.test$total.t + psi.3.hat*data.test$V.i*data.test$total.t)
  
  #Regret
  regret <- (Aopt.dwols - data.test$A.ijk)*(psi.0.hat + psi.1.hat*data.test$X.ijk + psi.2.hat*data.test$total.t + psi.3.hat*data.test$V.i*data.test$total.t)

  #E[Y|A = Aopt.dwols]
  #Optimal treatment based on what we observe in given dataset
  Aopt.dwols.obs <- as.numeric((psi.0.hat + psi.1.hat*data.sim$X.ijk + psi.2.hat*data.sim$total.t + psi.3.hat*data.sim$V.i*data.sim$total.t) > 0)
  
  #Mean of outcomes where treatment we observe is the same as the dwols optimal
  #treatment for our given dataset
  Y.agree <- mean(data.sim$Y.ijk[data.sim$A.ijk == Aopt.dwols.obs])
  
  return(data.frame(mean(Aopt.dwols == data.test$AOpt.ijk), mean(Y.opt), mean(Y.agree), psi.0.hat, psi.1.hat, psi.2.hat, psi.3.hat, xsi.0.hat, xsi.1.hat, mean(regret)))
}

#METHOD #4: ACCOUNTS FOR WITHIN- AND BETWEEN-GROUP INTERFERENCE, USES RESCALED NETWORK WEIGHTS FROM ALL NEIGHBOURS

dwols_mod_3 <- function(data.sim, data.test, I, J, n) {
  N <- I*J*n

  mod <- lm(Y.ijk ~ X.ijk + V.i + total.t + V.i:total.t + A.ijk + A.ijk:X.ijk + A.ijk:total.t + A.ijk:V.i:total.t, weights = w.L3.resc, data = data.sim)
  
  #Obtain extra blip parameter estimates from outcome model
  #(only coefficients for treatment terms)
  psi.0.hat <- coef(mod)[5]
  psi.1.hat <- coef(mod)[7]
  psi.2.hat <- coef(mod)[8]
  psi.3.hat <- coef(mod)[9]

  #Obtain spillover parameter estimates in the treatment-free model
  xsi.0.hat <- coef(mod)[4]
  xsi.1.hat <- coef(mod)[6]
  
  #Compute optimal treatment in testing set using blip parameter estimates
  Aopt.dwols <- as.numeric((psi.0.hat + psi.1.hat*data.test$X.ijk + psi.2.hat*data.test$total.t + psi.3.hat*data.test$V.i*data.test$total.t) > 0)
  
  #Optimal outcome
  Y.opt <- data.test$Y.ijk + (Aopt.dwols - data.test$A.ijk)*(psi.0.hat + psi.1.hat*data.test$X.ijk + psi.2.hat*data.test$total.t + psi.3.hat*data.test$V.i*data.test$total.t)
  
  #Regret
  regret <- (Aopt.dwols - data.test$A.ijk)*(psi.0.hat + psi.1.hat*data.test$X.ijk + psi.2.hat*data.test$total.t + psi.3.hat*data.test$V.i*data.test$total.t)

  #E[Y|A = Aopt.dwols]
  #Optimal treatment based on what we observe in given dataset
  Aopt.dwols.obs <- as.numeric((psi.0.hat + psi.1.hat*data.sim$X.ijk + psi.2.hat*data.sim$total.t + psi.3.hat*data.sim$V.i*data.sim$total.t) > 0)
  
  #Mean of outcomes where treatment we observe is the same as the dwols optimal
  #treatment for our given dataset
  Y.agree <- mean(data.sim$Y.ijk[data.sim$A.ijk == Aopt.dwols.obs])
  
  return(data.frame(mean(Aopt.dwols == data.test$AOpt.ijk), mean(Y.opt), mean(Y.agree), psi.0.hat, psi.1.hat, psi.2.hat, psi.3.hat, xsi.0.hat, xsi.1.hat, mean(regret)))
}

#******************************************************************************************************************************
#SIMULATIONS START HERE
Sims <- function(I, J, n, sam.size, tau3, tau2, sigma2, alpha.0, alpha.1, alpha.2, psi.0, psi.1, psi.2, psi.3, run, model) {
  #Storage for results
  Res <- data.frame(OPT = double(), DWOL = double(), DWOLI = double(), DWOLII = double(), DWOLIII = double(), YOPT = double(), YDWOL = double(), YDWOLI = double(), YDWOLII = double(), YDWOLIII = double(), YADWOL = double(), YADWOLI = double(), YADWOLII = double(), YADWOLIII = double(), REGRETDWOL = double(), REGRETDWOLI = double(), REGRETDWOLII = double(), REGRETDWOLIII = double())

  Min.Num.Alters <- c()
  Max.Num.Alters <- c()
  Mean.Num.Alters <- c()
  Median.Num.Alters <- c()

  Psi.0.Dwols <- c()

  Psi.1.Dwols <- c()

  Psi.2.Dwols <- c()

  Psi.3.Dwols <- c()

  Xsi.0.Dwols <- c()

  Xsi.1.Dwols <- c()


  Psi.0.Dwols.Mod.1 <- c()

  Psi.1.Dwols.Mod.1 <- c()

  Psi.2.Dwols.Mod.1 <- c()

  Psi.3.Dwols.Mod.1 <- c()

  Xsi.0.Dwols.Mod.1 <- c()

  Xsi.1.Dwols.Mod.1 <- c()


  Psi.0.Dwols.Mod.2 <- c()

  Psi.1.Dwols.Mod.2 <- c()

  Psi.2.Dwols.Mod.2 <- c()

  Psi.3.Dwols.Mod.2 <- c()

  Xsi.0.Dwols.Mod.2 <- c()

  Xsi.1.Dwols.Mod.2 <- c()


  Psi.0.Dwols.Mod.3 <- c()

  Psi.1.Dwols.Mod.3 <- c()

  Psi.2.Dwols.Mod.3 <- c()

  Psi.3.Dwols.Mod.3 <- c()

  Xsi.0.Dwols.Mod.3 <- c()

  Xsi.1.Dwols.Mod.3 <- c()



  for(seed in 1:run) {
    print(seed)
    #Generate data
    data.sim <- main.fun(I, J, n, sam.size, tau3, tau2, sigma2, alpha.0, alpha.1, alpha.2, psi.0, psi.1, psi.2, psi.3, seed)

    #Generate testing dataset
    data.test <- main.fun(I, J, n, sam.size, tau3, tau2, sigma2, alpha.0, alpha.1, alpha.2, psi.0, psi.1, psi.2, psi.3, seed+1)

    print(c("done"))

    #Extract summary information on number of alters in the test data
    Min.Num.Alters <- c(Min.Num.Alters, min(data.test$num.alters))
    Max.Num.Alters <- c(Max.Num.Alters, max(data.test$num.alters))
    Mean.Num.Alters <- c(Mean.Num.Alters, mean(data.test$num.alters))
    Median.Num.Alters <- c(Median.Num.Alters, median(data.test$num.alters))

    #Run each analysis
    res.dwols <- dwols(data.sim, data.test, I, J, n)
    res.dwols.mod.1 <- dwols_mod_1(data.sim, data.test, I, J, n)
    res.dwols.mod.2 <- dwols_mod_2(data.sim, data.test, I, J, n)
    res.dwols.mod.3 <- dwols_mod_3(data.sim, data.test, I, J, n)
    
    #Extract results from each method
    Aopt.dwols <- res.dwols[,1]
    Yopt.dwols <- res.dwols[,2]
    Yagree.dwols <- res.dwols[,3]
    Regret.dwols <- res.dwols[,10]

    psi.0.dwols <- res.dwols[,4]
    psi.1.dwols <- res.dwols[,5]
    psi.2.dwols <- res.dwols[,6]
    psi.3.dwols <- res.dwols[,7]

    xsi.0.dwols <- res.dwols[,8]
    xsi.1.dwols <- res.dwols[,9]

    Aopt.dwols.mod.1 <- res.dwols.mod.1[,1]
    Yopt.dwols.mod.1 <- res.dwols.mod.1[,2]
    Yagree.dwols.mod.1 <- res.dwols.mod.1[,3]
    Regret.dwols.mod.1 <- res.dwols.mod.1[,10]

    psi.0.dwols.mod.1 <- res.dwols.mod.1[,4]
    psi.1.dwols.mod.1 <- res.dwols.mod.1[,5]
    psi.2.dwols.mod.1 <- res.dwols.mod.1[,6]
    psi.3.dwols.mod.1 <- res.dwols.mod.1[,7]

    xsi.0.dwols.mod.1 <- res.dwols.mod.1[,8]
    xsi.1.dwols.mod.1 <- res.dwols.mod.1[,9]

    
    Aopt.dwols.mod.2 <- res.dwols.mod.2[,1]
    Yopt.dwols.mod.2 <- res.dwols.mod.2[,2]
    Yagree.dwols.mod.2 <- res.dwols.mod.2[,3]
    Regret.dwols.mod.2 <- res.dwols.mod.2[,10]

    psi.0.dwols.mod.2 <- res.dwols.mod.2[,4]
    psi.1.dwols.mod.2 <- res.dwols.mod.2[,5]
    psi.2.dwols.mod.2 <- res.dwols.mod.2[,6]
    psi.3.dwols.mod.2 <- res.dwols.mod.2[,7]

    xsi.0.dwols.mod.2 <- res.dwols.mod.2[,8]
    xsi.1.dwols.mod.2 <- res.dwols.mod.2[,9]
    
    Aopt.dwols.mod.3 <- res.dwols.mod.3[,1]
    Yopt.dwols.mod.3 <- res.dwols.mod.3[,2]
    Yagree.dwols.mod.3 <- res.dwols.mod.3[,3]
    Regret.dwols.mod.3 <- res.dwols.mod.3[,10]

    psi.0.dwols.mod.3 <- res.dwols.mod.3[,4]
    psi.1.dwols.mod.3 <- res.dwols.mod.3[,5]
    psi.2.dwols.mod.3 <- res.dwols.mod.3[,6]
    psi.3.dwols.mod.3 <- res.dwols.mod.3[,7]

    xsi.0.dwols.mod.3 <- res.dwols.mod.3[,8]
    xsi.1.dwols.mod.3 <- res.dwols.mod.3[,9]

    
    #True optimal outcome
    Yopt <- data.test$Y.ijk + (data.test$AOpt.ijk - data.test$A.ijk)*(psi.0 + psi.1*data.test$X.ijk + psi.2*data.test$total.t + psi.3*data.test$V.i*data.test$total.t)
    
    #Store results
    Res[seed,] <- c(mean(data.sim$A.ijk == data.sim$AOpt.ijk), Aopt.dwols, Aopt.dwols.mod.1, Aopt.dwols.mod.2, Aopt.dwols.mod.3, mean(Yopt), Yopt.dwols, Yopt.dwols.mod.1, Yopt.dwols.mod.2, Yopt.dwols.mod.3, Yagree.dwols, Yagree.dwols.mod.1, Yagree.dwols.mod.2, Yagree.dwols.mod.3, Regret.dwols, Regret.dwols.mod.1, Regret.dwols.mod.2, Regret.dwols.mod.3)
    Psi.0.Dwols <- c(Psi.0.Dwols, psi.0.dwols)
    Psi.1.Dwols <- c(Psi.1.Dwols, psi.1.dwols)
    Psi.2.Dwols <- c(Psi.2.Dwols, psi.2.dwols)
    Psi.3.Dwols <- c(Psi.3.Dwols, psi.3.dwols)

    Xsi.0.Dwols <- c(Xsi.0.Dwols, xsi.0.dwols)
    Xsi.1.Dwols <- c(Xsi.1.Dwols, xsi.1.dwols)

    Psi.0.Dwols.Mod.1 <- c(Psi.0.Dwols.Mod.1, psi.0.dwols.mod.1)
    Psi.1.Dwols.Mod.1 <- c(Psi.1.Dwols.Mod.1, psi.1.dwols.mod.1)
    Psi.2.Dwols.Mod.1 <- c(Psi.2.Dwols.Mod.1, psi.2.dwols.mod.1)
    Psi.3.Dwols.Mod.1 <- c(Psi.3.Dwols.Mod.1, psi.3.dwols.mod.1)

    Xsi.0.Dwols.Mod.1 <- c(Xsi.0.Dwols.Mod.1, xsi.0.dwols.mod.1)
    Xsi.1.Dwols.Mod.1 <- c(Xsi.1.Dwols.Mod.1, xsi.1.dwols.mod.1)


    Psi.0.Dwols.Mod.2 <- c(Psi.0.Dwols.Mod.2, psi.0.dwols.mod.2)
    Psi.1.Dwols.Mod.2 <- c(Psi.1.Dwols.Mod.2, psi.1.dwols.mod.2)
    Psi.2.Dwols.Mod.2 <- c(Psi.2.Dwols.Mod.2, psi.2.dwols.mod.2)
    Psi.3.Dwols.Mod.2 <- c(Psi.3.Dwols.Mod.2, psi.3.dwols.mod.2)

    Xsi.0.Dwols.Mod.2 <- c(Xsi.0.Dwols.Mod.2, xsi.0.dwols.mod.2)
    Xsi.1.Dwols.Mod.2 <- c(Xsi.1.Dwols.Mod.2, xsi.1.dwols.mod.2)

    Psi.0.Dwols.Mod.3 <- c(Psi.0.Dwols.Mod.3, psi.0.dwols.mod.3)
    Psi.1.Dwols.Mod.3 <- c(Psi.1.Dwols.Mod.3, psi.1.dwols.mod.3)
    Psi.2.Dwols.Mod.3 <- c(Psi.2.Dwols.Mod.3, psi.2.dwols.mod.3)
    Psi.3.Dwols.Mod.3 <- c(Psi.3.Dwols.Mod.3, psi.3.dwols.mod.3)

    Xsi.0.Dwols.Mod.3 <- c(Xsi.0.Dwols.Mod.3, xsi.0.dwols.mod.3)
    Xsi.1.Dwols.Mod.3 <- c(Xsi.1.Dwols.Mod.3, xsi.1.dwols.mod.3)


  }

  Res.List <- list("Res" = Res, "Psi.0.Dwols" = Psi.0.Dwols, "Psi.1.Dwols" = Psi.1.Dwols, "Psi.2.Dwols" = Psi.2.Dwols, "Psi.3.Dwols" = Psi.3.Dwols, "Xsi.0.Dwols" = Xsi.0.Dwols, "Xsi.1.Dwols" = Xsi.1.Dwols,"Psi.0.Dwols.Mod.1" = Psi.0.Dwols.Mod.1, "Psi.1.Dwols.Mod.1" = Psi.1.Dwols.Mod.1, "Psi.2.Dwols.Mod.1" = Psi.2.Dwols.Mod.1, "Psi.3.Dwols.Mod.1" = Psi.3.Dwols.Mod.1, "Xsi.0.Dwols.Mod.1" = Xsi.0.Dwols.Mod.1, "Xsi.1.Dwols.Mod.1" = Xsi.1.Dwols.Mod.1, "Psi.0.Dwols.Mod.2" = Psi.0.Dwols.Mod.2, "Psi.1.Dwols.Mod.2" = Psi.1.Dwols.Mod.2, "Psi.2.Dwols.Mod.2" = Psi.2.Dwols.Mod.2, "Psi.3.Dwols.Mod.2" = Psi.3.Dwols.Mod.2, "Xsi.0.Dwols.Mod.2" = Xsi.0.Dwols.Mod.2, "Xsi.1.Dwols.Mod.2" = Xsi.1.Dwols.Mod.2,"Psi.0.Dwols.Mod.3" = Psi.0.Dwols.Mod.3, "Psi.1.Dwols.Mod.3" = Psi.1.Dwols.Mod.3, "Psi.2.Dwols.Mod.3" = Psi.2.Dwols.Mod.3, "Psi.3.Dwols.Mod.3" = Psi.3.Dwols.Mod.3, "Xsi.0.Dwols.Mod.3" = Xsi.0.Dwols.Mod.3, "Xsi.1.Dwols.Mod.3" = Xsi.1.Dwols.Mod.3, "Min.Num.Alters" = Min.Num.Alters, "Max.Num.Alters" = Max.Num.Alters, "Mean.Num.Alters" = Mean.Num.Alters, "Median.Num.Alters" = Median.Num.Alters)
  return(Res.List)
}

#Perform "run" number of simulation runs
#Sample size of "n" patients per Level 2 group
#Sample size of egos is denoted "sam.size"
#Vary the no. of Level 2 groups: "a" lowest, "b" highest, sequence by 1
#c refers to number of Level 3 groups
begin.Sims <- function(run, n, sam.size, a, b, c, tau3, tau2, sigma2, alpha.0, alpha.1, alpha.2, psi.0, psi.1, psi.2, psi.3){
  
  J.list <- seq(a, b, 1)
  
  #Storage for plotting later
  A.Opts <- A.Opts.D <- A.Opts.1 <- A.Opts.2 <- A.Opts.3 <-  matrix(nrow = run, ncol = length(J.list))
  Y.Opts <- Y.Opts.D <- Y.Opts.1 <- Y.Opts.2 <- Y.Opts.3 <-  matrix(nrow = run, ncol = length(J.list))
  Y.Agree <- Y.Agree.D <- Y.Agree.1 <- Y.Agree.2 <- Y.Agree.3 <-  matrix(nrow = run, ncol = length(J.list))
  Regret.D <- Regret.1 <- Regret.2 <- Regret.3 <-  matrix(nrow = run, ncol = length(J.list))

  Min.Alters <- matrix(nrow=run, ncol=length(J.list))

  Max.Alters <- matrix(nrow=run, ncol=length(J.list))

  Mean.Alters <- matrix(nrow=run, ncol=length(J.list))

  Median.Alters <- matrix(nrow=run, ncol=length(J.list))


  Psi.0.Dw  <- matrix(nrow = run, ncol = length(J.list))

  Psi.1.Dw  <- matrix(nrow = run, ncol = length(J.list))

  Psi.2.Dw  <- matrix(nrow = run, ncol = length(J.list))

  Psi.3.Dw  <- matrix(nrow = run, ncol = length(J.list))
  
  Xsi.0.Dw <- matrix(nrow=run, ncol=length(J.list))

  Xsi.1.Dw <- matrix(nrow=run, ncol=length(J.list))


  Psi.0.Dw.M1  <- matrix(nrow = run, ncol = length(J.list))

  Psi.1.Dw.M1  <- matrix(nrow = run, ncol = length(J.list))

  Psi.2.Dw.M1  <- matrix(nrow = run, ncol = length(J.list))

  Psi.3.Dw.M1  <- matrix(nrow = run, ncol = length(J.list))

  Xsi.0.Dw.M1  <- matrix(nrow = run, ncol = length(J.list))

  Xsi.1.Dw.M1  <- matrix(nrow = run, ncol = length(J.list))


  Psi.0.Dw.M2  <- matrix(nrow = run, ncol = length(J.list))

  Psi.1.Dw.M2  <- matrix(nrow = run, ncol = length(J.list))

  Psi.2.Dw.M2  <- matrix(nrow = run, ncol = length(J.list))

  Psi.3.Dw.M2  <- matrix(nrow = run, ncol = length(J.list))

  Xsi.0.Dw.M2  <- matrix(nrow = run, ncol = length(J.list))

  Xsi.1.Dw.M2  <- matrix(nrow = run, ncol = length(J.list))


  Psi.0.Dw.M3  <- matrix(nrow = run, ncol = length(J.list))

  Psi.1.Dw.M3  <- matrix(nrow = run, ncol = length(J.list))

  Psi.2.Dw.M3  <- matrix(nrow = run, ncol = length(J.list))

  Psi.3.Dw.M3  <- matrix(nrow = run, ncol = length(J.list))

  Xsi.0.Dw.M3  <- matrix(nrow = run, ncol = length(J.list))

  Xsi.1.Dw.M3  <- matrix(nrow = run, ncol = length(J.list))
  
  #Set up counter
  count <- 0
  
  for (j in J.list) {
    count <- count + 1
    print(j)
    #Run simulations
    Res.List <- Sims(c, j, n, sam.size, tau3, tau2, sigma2, alpha.0, alpha.1, alpha.2, psi.0, psi.1, psi.2, psi.3, run, model)
    Res <- Res.List$Res

    #Store optimal treatment rate for each method
    A.Opts[,count] <- Res[,1]
    A.Opts.D[,count] <- Res[,2]
    A.Opts.1[,count] <- Res[,3]
    A.Opts.2[,count] <- Res[,4]
    A.Opts.3[,count] <- Res[,5]    
    
    #Store mean optimal outcome
    Y.Opts[,count] <- Res[,6]
    Y.Opts.D[,count] <- Res[,7]
    Y.Opts.1[,count] <- Res[,8]
    Y.Opts.2[,count] <- Res[,9]
    Y.Opts.3[,count] <- Res[,10]    
  
    #Store "agreements"
    Y.Agree.D[,count] <- Res[,11]
    Y.Agree.1[,count] <- Res[,12]
    Y.Agree.2[,count] <- Res[,13]
    Y.Agree.3[,count] <- Res[,14]  

    #Store regrets
    Regret.D[,count] <- Res[,15]
    Regret.1[,count] <- Res[,16]
    Regret.2[,count] <- Res[,17]
    Regret.3[,count] <- Res[,18] 

    #Store min, max, mean, and median alters in each simulation run

    Min.Alters[,count] <- Res.List$Min.Num.Alters  
    Max.Alters[,count] <- Res.List$Max.Num.Alters
    Mean.Alters[,count] <- Res.List$Mean.Num.Alters
    Median.Alters[,count] <- Res.List$Median.Num.Alters

    #Store Psi Estimates under each method

    Psi.0.Dw[,count]  <- Res.List$Psi.0.Dwols

    Psi.1.Dw[,count]  <- Res.List$Psi.1.Dwols

    Psi.2.Dw[,count]  <- Res.List$Psi.2.Dwols

    Psi.3.Dw[,count]  <- Res.List$Psi.3.Dwols

    Xsi.0.Dw[,count]  <- Res.List$Xsi.0.Dwols

    Xsi.1.Dw[,count]  <- Res.List$Xsi.1.Dwols



    Psi.0.Dw.M1[,count]  <- Res.List$Psi.0.Dwols.Mod.1

    Psi.1.Dw.M1[,count]  <- Res.List$Psi.1.Dwols.Mod.1

    Psi.2.Dw.M1[,count]  <- Res.List$Psi.2.Dwols.Mod.1

    Psi.3.Dw.M1[,count]  <- Res.List$Psi.3.Dwols.Mod.1

    Xsi.0.Dw.M1[,count]  <- Res.List$Xsi.0.Dwols.Mod.1

    Xsi.1.Dw.M1[,count]  <- Res.List$Xsi.1.Dwols.Mod.1



    Psi.0.Dw.M2[,count]  <- Res.List$Psi.0.Dwols.Mod.2

    Psi.1.Dw.M2[,count]  <- Res.List$Psi.1.Dwols.Mod.2

    Psi.2.Dw.M2[,count]  <- Res.List$Psi.2.Dwols.Mod.2

    Psi.3.Dw.M2[,count]  <- Res.List$Psi.3.Dwols.Mod.2

    Xsi.0.Dw.M2[,count]  <- Res.List$Xsi.0.Dwols.Mod.2

    Xsi.1.Dw.M2[,count]  <- Res.List$Xsi.1.Dwols.Mod.2


    Psi.0.Dw.M3[,count]  <- Res.List$Psi.0.Dwols.Mod.3

    Psi.1.Dw.M3[,count]  <- Res.List$Psi.1.Dwols.Mod.3

    Psi.2.Dw.M3[,count]  <- Res.List$Psi.2.Dwols.Mod.3

    Psi.3.Dw.M3[,count]  <- Res.List$Psi.3.Dwols.Mod.3

    Xsi.0.Dw.M3[,count]  <- Res.List$Xsi.0.Dwols.Mod.3

    Xsi.1.Dw.M3[,count]  <- Res.List$Xsi.1.Dwols.Mod.3

    
    #Use counter to mark progress
    cat('count ', count, ' done\n')
  }

  #Stack values from each method
  #Each of these variables: rows=run*5, cols = length(J.list)
  A<- c(rbind(A.Opts, A.Opts.D, A.Opts.1, A.Opts.2, A.Opts.3))
  YO <- c(rbind(Y.Opts, Y.Opts.D, Y.Opts.1, Y.Opts.2, Y.Opts.3))
  YA <- c(rbind(Y.Agree, Y.Agree.D, Y.Agree.1, Y.Agree.2, Y.Agree.3))
  #rows=run*4, cols = length(J.list)
  RegretResults <- c(rbind(Regret.D, Regret.1, Regret.2, Regret.3))

  print(A)
  print(YO)
  print(YA)
  print(RegretResults)


  mthd <- rep(c("Observed", "ME1", "ME2", "ME3", "ME4"), each=run)
  mthd_regret <- rep(c("ME1", "ME2", "ME3", "ME4"), each=run)
  psi.type.dw <- rep(c("psi.0", "psi.1", "psi.2", "psi.3"), each=run)
  psi.type.1 <- rep(c("psi.0", "psi.1", "psi.2", "psi.3"), each=run)
  psi.type.2 <- rep(c("psi.0", "psi.1", "psi.2", "psi.3"), each=run)
  psi.type.3 <- rep(c("psi.0", "psi.1", "psi.2", "psi.3"), each=run)

  xsi.type.dw <- rep(c("xsi.0", "xsi.1"), each=run)
  xsi.type.1 <- rep(c("xsi.0", "xsi.1"), each=run)
  xsi.type.2 <- rep(c("xsi.0", "xsi.1"), each=run)
  xsi.type.3 <- rep(c("xsi.0", "xsi.1"), each=run)


  AOpts <- data.frame(mthd, as.numeric(A))
  names(AOpts) <- c("Method", "TrtRate")
  YOpts <- data.frame(mthd, as.numeric(YO))
  names(YOpts) <- c("Method", "Outcome")
  YAgree <- data.frame(mthd, as.numeric(YA))
  names(YAgree) <- c("Method", "Agreement")
  RegretMethods <- data.frame(mthd_regret, as.numeric(RegretResults))
  names(RegretMethods) <- c("Method", "Regret")

  Min.Alters <- as.data.frame(Min.Alters)
  names(Min.Alters) <- as.character(J.list)
  Max.Alters <- as.data.frame(Max.Alters)
  names(Max.Alters) <- as.character(J.list)
  Mean.Alters <- as.data.frame(Mean.Alters)
  names(Mean.Alters) <- as.character(J.list)
  Median.Alters <- as.data.frame(Median.Alters)
  names(Median.Alters) <- as.character(J.list)

  Psi.Dw <- data.frame(psi.type.dw, as.numeric(c(rbind(Psi.0.Dw, Psi.1.Dw, Psi.2.Dw, Psi.3.Dw))))
  names(Psi.Dw) <- c("PsiParameter", "Estimate")

  Xsi.Dw <- data.frame(xsi.type.dw, as.numeric(c(rbind(Xsi.0.Dw, Xsi.1.Dw))))
  names(Xsi.Dw) <- c("XsiParameter", "Estimate")

  Psi.M1 <- data.frame(psi.type.1, as.numeric(c(rbind(Psi.0.Dw.M1, Psi.1.Dw.M1, Psi.2.Dw.M1, Psi.3.Dw.M1))))
  names(Psi.M1) <- c("PsiParameter", "Estimate")

  Xsi.M1 <- data.frame(xsi.type.1, as.numeric(c(rbind(Xsi.0.Dw.M1, Xsi.1.Dw.M1))))
  names(Xsi.M1) <- c("XsiParameter", "Estimate")

  Psi.M2 <- data.frame(psi.type.2, as.numeric(c(rbind(Psi.0.Dw.M2, Psi.1.Dw.M2, Psi.2.Dw.M2, Psi.3.Dw.M2))))
  names(Psi.M2) <- c("PsiParameter", "Estimate")

  Xsi.M2 <- data.frame(xsi.type.2, as.numeric(c(rbind(Xsi.0.Dw.M2, Xsi.1.Dw.M2))))
  names(Xsi.M2) <- c("XsiParameter", "Estimate")

  Psi.M3 <- data.frame(psi.type.3, as.numeric(c(rbind(Psi.0.Dw.M3, Psi.1.Dw.M3, Psi.2.Dw.M3, Psi.3.Dw.M3))))
  names(Psi.M3) <- c("PsiParameter", "Estimate")

  Xsi.M3 <- data.frame(xsi.type.3, as.numeric(c(rbind(Xsi.0.Dw.M3, Xsi.1.Dw.M3))))
  names(Xsi.M3) <- c("XsiParameter", "Estimate")


  bp_list <- list("AOpts" = AOpts, "YOpts" = YOpts, "YAgree" = YAgree, "RegretMethods" = RegretMethods, "Psi.Dw" = Psi.Dw, "Psi.M1" = Psi.M1, "Psi.M2" = Psi.M2, "Psi.M3" = Psi.M3, "Xsi.Dw" = Xsi.Dw, "Xsi.M1" = Xsi.M1, "Xsi.M2" = Xsi.M2, "Xsi.M3" = Xsi.M3, "Min.Alters" = Min.Alters, "Max.Alters" = Max.Alters, "Mean.Alters" = Mean.Alters, "Median.Alters" = Median.Alters)
  #bp_list <- list("AOpts" = AOpts, "YOpts" = YOpts, "YAgree" = YAgree)

  
  return(bp_list)
}

#***************************************************************************************************
