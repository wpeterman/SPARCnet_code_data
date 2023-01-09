## Code to conduct population projection simulation
## Bill Peterman
## 23 June 2020
## Updated 9 January 2023
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Libraries ---------------------------------------------------------------

library(jagsUI)
library(brms)

# Data --------------------------------------------------------------------

## Import growth model results
## Posterior parameter estimates from fitted growth model
growth_df <- readRDS("fitted_model_df.rds")

## Import Clutch size results from 
## Posterior parameter estimates from linear regression of Lotter data
egg_df <- read.csv("egg_model_normal_posteriors.csv")[1:nrow(growth_df),]


## Data used to fit regression between SVL and clutch size
# lotter <- read.csv("Lotter_Data.csv")


# Simulation --------------------------------------------------------------

n_ind <- 100000 # Number of individuals
yrs <- 20 + 1

## Generate normally distributed survival values based on Munoz
## Estimates from our spatial capture-recap model unrealistically high
phi <- TruncatedNormal::rtnorm(1000, 0.836, sd = 0.07, lb = 0.4, ub = 1)
hist(phi)

off <- svl <- surv <- array(NA, c(n_ind, yrs))
out_list <- vector('list', 4)
names(out_list) <- c('slope', 'ridge',
                     'slope_m', 'ridge_m')

for(i in 1:4){
  if(i == 1){ # Slope female
    sex <- 2
    pos_sex <- 7
    pos <- 3
    male <- 0
  } else if(i == 2)  { # Ridge female
    sex <- 2
    pos_sex <- 9
    pos <- 4
    male <- 0
  } else if(i == 3){ # Slope Male
    sex <- 1
    pos_sex <- 6 
    pos <- 3
    male <- 1
  } else {       # Ridge Male
    sex <- 1
    pos_sex <- 8
    pos <- 4
    male <- 1
  }
  
  for(y in 1:yrs){
    samp <- sort(sample(1:nrow(egg_df), n_ind, replace = F))
    
    if(y == 1){
      surv[,y] <- 1
      svl[,y] <- 13.5
      off[,y] <- 0
    } else {
      phi <- TruncatedNormal::rtnorm(n_ind, 0.836, sd = 0.07, lb = 0.5, ub = 1)
      surv[,y] <- surv[,y-1] * rbinom(n_ind, 1, phi)
      
      ## Apply growth model
      K_ <- exp(growth_df[samp,pos] + growth_df$K_male[samp] * male)
      ## WORKS
      svl[,y] <- rnorm(n = n_ind, svl[,y-1], 1/growth_df$sigma.SVL[samp]^2) + (growth_df[samp,sex] - svl[,y-1]) * (1 - exp(-K_))
      
      ## Determine if individual is gravid
      ## If so, how many eggs
      gravid <- rep(0, length(svl[,1]))
      gravid[svl[,y] >= 34 & svl[,y] < 43] <- rbinom(sum(svl[,y] >= 34 & svl[,y] < 43), 1, prob = 0.56)
      gravid[svl[,y] >= 43] <- rbinom(sum(svl[,y] >= 43), 1, prob = 0.94)
      
      ## Of eggs laid, how many survive
      off[,y] <- floor(egg_df[samp,1] + egg_df[samp, 2] * svl[,y] * 0.9) * gravid * surv[,y]
      
    } 
  } # End year loop
  out_list[[i]] <- list(surv = surv,
                        svl = svl,
                        off = off)
} # End location loop


# Summarize Results -------------------------------------------------------
# * Age -------------------------------------------------------------------

# >> Slope ----------------------------------------------------------------
age_s <- apply(out_list$slope$surv, 1, function(x) which(x == 0)[1] - 1)
age_s[is.na(age_s)] <- yrs - 1
hist(age_s)
(med_age_s <- median(age_s))
(mean_age_s <- mean(age_s))
(ci_age_s <- quantile(age_s, c(0.025, 0.975)))
(sd_age_s <- sd(age_s))

# >> Ridge ----------------------------------------------------------------
age_r <- apply(out_list$ridge$surv, 1, function(x) which(x == 0)[1] - 1)
age_r[is.na(age_r)] <- yrs - 1
hist(age_r)
(med_age_r <- median(age_r))
(mean_age_r <- mean(age_r))
(ci_age_r <- quantile(age_r, c(0.025, 0.975)))
(sd_age_r <- sd(age_r))


# * Clutches --------------------------------------------------------------
# >> Slope ----------------------------------------------------------------
n_clutches_s <- apply(out_list$slope$off, 1, function(x) sum(x != 0))
hist(n_clutches_s)
(med_clutch_s <- median(n_clutches_s))
(mn_clutch_s <- mean(n_clutches_s))
(sd_clutch_s <- sd(n_clutches_s))
(ci_clutch_s <- quantile(n_clutches_s, c(0.025, 0.975)))

# >> Ridge ----------------------------------------------------------------
n_clutches_r <- apply(out_list$ridge$off, 1, function(x) sum(x != 0))
hist(n_clutches_r)
(med_clutch_r <- median(n_clutches_r))
(mn_clutch_r <- mean(n_clutches_r))
(sd_clutch_r <- sd(n_clutches_r))
(ci_clutch_r <- quantile(n_clutches_r, c(0.025, 0.975)))


# * Number of Offspring ---------------------------------------------------
# >> Slope ----------------------------------------------------------------
t_off_s <- apply(out_list$slope$off, 1, sum)
hist(t_off_s)
(med_off_s <- median(t_off_s))
(mn_off_s <- mean(t_off_s))
(sd_off_s <- sd(t_off_s))
(hmn_off_s <- psych::harmonic.mean(t_off_s, zero = F))
(ci_off_s <- quantile(t_off_s, c(0.025, 0.975)))

# >> Ridge ----------------------------------------------------------------
t_off_r <- apply(out_list$ridge$off, 1, sum)
hist(t_off_r)
(med_off_r <- median(t_off_r))
(mn_off_r <- mean(t_off_r))
(sd_off_r <- sd(t_off_r))
(hmn_off_r <- psych::harmonic.mean(t_off_r, zero = F))
(ci_off_r <- quantile(t_off_r, c(0.025, 0.975)))


# Fecundity increase ------------------------------------------------------

((mn_off_r - mn_off_s) / mn_off_s) * 100
((hmn_off_r - hmn_off_s) / hmn_off_s) * 100