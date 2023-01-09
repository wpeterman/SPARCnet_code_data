### Info ######################################################################################
## Code to fit growth model to mark-recapture data
## Written by Bill Peterman
## Updated 9 January 2023
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Prep data ---------------------------------------------------------------


## Load Libraries
library(jagsUI)
library(readr)

## Import SVL data
# spreadsheet with: date, survey_number, R_S, plot, salamander ID, sex, SVL
DAT <- read.csv("data/SalDatGrowth.csv", header = TRUE)


## Create Individual x survey array of SVLs (Y)
NSurv <- max(DAT$SurveyNum)
NInd <- length(unique(DAT$UniqueID))
Y <- matrix(nrow=NInd,ncol=NSurv)  # create empty capture history matrix
Plot <- numeric()

### Capture histories -------------------------------------------------------

iter <- 0   # 'iter' will be an object to track which row of the capture history to fill

for (i in unique(DAT$UniqueID)){   # cycle through the codes
  iter <- iter+1                       # add 1 to our capture history row number tracker
  caps <- which(DAT$UniqueID==i)   # which rows in database have the current code
  Plot[iter] <- DAT$PosID[caps[1]]   # extract plot at first capture
  for (j in caps){              # loop over rows with capture record records for current individual
    surv <- DAT$SurveyNum[j]
    svl <- DAT$SVL[j]
    Y[iter,surv] <- svl
  }
}

## Create vector of first capture occasions (f)
get.first <- function(x) min(which(x!=0))
f <- apply(Y,1,get.first)

f <- f[which(is.finite(f))]  # remove individuals with no SVLs measured - [none here]
Y <- Y[which(is.finite(f)),] # and remove them here too

## Import matrix of SURVEY DATES and convert to INTERVALS between sequential surveys
Dates2 <- t(read.csv("data/DateMat_osu.csv", header=T)[,-c(1,2)])


SampleDates <- julian(as.Date(Dates2, format="%m/%d/%Y"))
SampleDates <- matrix(SampleDates,nrow=4,byrow=T) 
# nrow is the number of rows in the spreadsheet (aka, number of plots)
SampleDates <- SampleDates[,-c(7,8)]
SampleDates <- t(apply(SampleDates,MARGIN=1,FUN=diff))

Interval <- cbind(rep(NA,4),SampleDates) # 4 = number of plots

## Create vector of last survey occasions for each plot (l)
get.last <- function(x) max(which(!is.na(x)))
l <- apply(Interval,1,get.last)
l <- l[Plot]


### Covariates --------------------------------------------------------------

## Sex
sex.df <- read.csv("data/SexMat.csv")
sex.mat <- as.matrix(sex.df[,-c(1)])

## RS [factor]
rs.df <- read.csv("data/RSmat.csv")
rs.mat <- as.matrix(rs.df[,-c(1)])

## RS [random effect] 
pos.mat <- rs.mat + 1 # 1 = slope, 2 = ridge


#### MODEL - RS + Sex --------------------------------------------------------------
# RS as a random effect

###~~ Write model -------------------------------------------------------------

sink("growth_mods/Growth_pos_sex_FULL.txt") # sink is output file
cat("
model{
  for (i in 1:NInd){     # Looping over all capture records in order (i.e., no need to format as a capture history)
    for (j in f[i]:f[i]){
      #Fabens mark-recapture formulation
      SVL[i,j] ~ dnorm(SVL0[i,j], tau.SVL)  # SVL[i] is measured size (i.e., your response)
      SVL0[i,j] ~ dunif(10,60)     # SVL0 is the true size at capture
    }
    for (j in (f[i]+1):l[i]){
      #Fabens mark-recapture formulation
      SVL[i,j] ~ dnorm(SVL0[i,j], tau.SVL)  # SVL[i] is measured size (i.e., your response)
     
      SVL0[i,j] <- SVL0[i,j-1] +  (L.inf[sex[i,j]] - SVL0[i,j-1]) * (1 - exp(-K[i,j]*(Interval[Plot[i],j]/365) ))
      
      log(K[i,j]) <- K_pos[pos[i,j]] + 
        K_m * (2 - sex[i,j]) 
    }
  }
  
  tau.SVL ~ dgamma(1,0.001)  # prior on precision
  sigma.SVL <- sqrt(1/tau.SVL)  # converting precision to sd
  
  # growth coefficient
  K_m ~ dnorm(0,0.01)
  K_pre ~ dnorm(0,0.01)
  L_slope ~ dnorm(0,0.01)
  
  # Prior Random site effect (Group by site)
  for(i in 1:2) {
    K_pos[i] ~ dnorm(0, 0.01)      # Position random effects
    K_sex[i] ~ dnorm(0, 0.01)      # Sex effect
    L.inf[i] ~ dnorm(48, 0.01)T(0,)   # Prior for asymptotic size parameter
  }
 
  ## Derived param
  K_r_m <- K_pos[2] + K_m 
  K_r_f <- K_pos[2]
  K_s_m <- K_pos[1] + K_m 
  K_s_f <- K_pos[1]
 
} # close model loop
",fill = TRUE)

sink()

## Bundle data
recap <- which(apply(sex.mat, 1, function(x) any(x != 0)))
sex_mat_final <- (sex.mat[recap,1:6]) + ((sex.mat[recap,7:12] + 1) * (sex.mat[recap,7:12]))
sum(apply(sex_mat_final, 1, function(x) any(x > 2)))
sum(apply(sex_mat_final, 1, function(x) any(x == 0)))
recap2 <- recap
recap2 <- recap[-which(apply(sex_mat_final, 1, function(x) any(x == 0) | any(x > 2)))]
sex_ <- sex_mat_final[-which(apply(sex_mat_final, 1, function(x) any(x == 0) | any(x > 2))),]
SVL <- Y[recap2,]; NInd <- nrow(SVL)

jags.data <- list(SVL = SVL,
                  NInd = NInd,
                  f = f[recap2],
                  l = l[recap2],
                  Plot = Plot[recap2],
                  Interval = Interval,
                  pos = pos.mat[recap2,],
                  pre = as.vector(scale(surv_dates$pre.sum)),
                  sex = sex_)

inits <- NULL

parameters <- c("L.inf", 
                "K_pos", 
                "K_pre", 
                'K_m',
                'K_r_m',
                "K_r_f",
                "K_s_m",
                'K_s_f',
                'sigma.SVL',
                'tau.SVL',
                'Sf_r.a',
                'Sf_r.b',
                'Sm_r.a',
                'Sm_r.b')


###~~ Run model ---------------------------------------------------------------

system.time(OUT <- jags(data=jags.data, inits=inits,
                        parameters, 
                        model.file="growth_mods/Growth_pos_sex_FULL.txt", parallel = T,
                        n.chains = 5,
                        n.adapt = 5000,
                        n.iter = 200000,
                        n.burnin = 25000,
                        n.thin = 5))

# saveRDS(OUT, file = "growth_mods/fitted_model.rds")


mod_df <- data.frame(OUT$sims.list)
names(mod_df) <- c('size_male',
                   'size_female',
                   'K_slope',
                   'K_ridge',
                   'K_male',
                   'K_ridge_m',
                   'K_ridge_f',
                   'K_slope_m',
                   'K_slope_f',
                   'deviance')
# saveRDS(mod_df, file = "ProjectionSim/fitted_model_df.rds")

mean(mod_df$size_male < mod_df$size_female)
mean(mod_df$K_slope < mod_df$K_ridge)

mean(exp(mod_df$K_slope)); mean(exp(mod_df$K_ridge))
mean(exp(mod_df$K_ridge_m)); mean(exp(mod_df$K_slope_m))
mean(exp(mod_df$K_ridge_f)); mean(exp(mod_df$K_slope_f))

mean(exp(mod_df$K_ridge_m) > exp(mod_df$K_slope_m)) # Males between positions
mean(exp(mod_df$K_ridge_f) > exp(mod_df$K_slope_f)) # Females between positions
mean(exp(mod_df$K_slope) < exp(mod_df$K_ridge)) # Position
mean(exp(mod_df$K_ridge_m) > exp(mod_df$K_ridge_f)) # Sex within ridge
mean(exp(mod_df$K_slope_m) > exp(mod_df$K_slope_f)) # Sex within slope