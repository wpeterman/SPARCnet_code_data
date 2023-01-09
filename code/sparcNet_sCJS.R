## Code to fit spatial capture-recapture model to 
## red-backed salamander mark-recap data
## 
## Written by Bill Peterman
## Updated 9 January 2023
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(jagsUI)
# Load Data ---------------------------------------------------------------

## Import data needed to fit model
data0 <- readRDS("data/sCJS_Data.rds")

# Jags Model (NULL) --------------------------------------------------------------
  
writeLines(
  "
  model{
    ## PRIORS AND CONSTRAINTS: ##
    #Space use and recapture probability parameters:
    for(s in 1:2){
      sigma[s] ~ dunif(0.1,20)
    }
    for(s in 1:2){ #Creates additive effect of being male on detection process
      lambda[s] <- exp(log.lambda0) * pow(beta,(s-1))
    }
    PL ~ dunif(0.01,0.99)
    log.lambda0 <- log(-log(1-PL))
    beta ~ dunif(0.1,10)
    
    
    ##---- NULL  MOD--
    for(s in 1:2){
      # Phi[s] ~ dlogis(0,1) #survival for each site
      Phi[s] ~ dbeta(9,3) #survival for each site (probability scale)
      
      # Phi[s] <- logit(pPhi[s]) # Convert to logit scale
      
      for(k in 1:(n.prim - 1)){
        phi[s,k] <- pow(Phi[s], dt[k])
      }
    }
    ##--- END  NULL MOD  
    
    # Dispersal parameters:
    for(s in 1:2){
      dmean[s] ~ dunif(0,100) #mean dispersal distance
      dlambda[s] <- 1/dmean[s]
    }
    ## MODEL: ##
    #Loop over individuals that are only seen in the last primary session or the
    #session they are censored
    for(i in 1:N[1]){
      z[i,first[i]] ~ dbern(1)
      S[i,1,first[i]] ~ dunif(xlow[i], xupp[i]) # Prior for the first x coordinate
      S[i,2,first[i]] ~ dunif(ylow[i], yupp[i]) # Prior for the first y coordinate
      g[i,first[i],1] <- 0
      for(r in 1:R){ # trap
        D[i,r,first[i]] <- sqrt(pow(S[i,1,first[i]]-X[r,1],2) + pow(S[i,2,first[i]]-X[r,2],2)) #Encounter model
        g[i,first[i],r+1] <- exp(-pow(D[i,r,first[i]]/sigma[gr[i]], 2)) # Trap exposure
      }
      G[i,first[i]] <- sum(g[i,first[i],]) # Total trap exposure
      for(j in 1:J[i,first[i]]){
        P[i,j,first[i]] <- 1 - exp(-lambda[gr[i]]*G[i,first[i]]) # Probability of being captured
        PI[i,first[i],j] <- step(H[i,j,first[i]]-2)*(g[i,first[i],H[i,j,first[i]]]/(G[i,first[i]]+ 0.0001))*P[i,j,first[i]] + (1.0001-step(H[i,j,first[i]]-2))*(1.0001-P[i,j,first[i]])
        Ones[i,j,first[i]] ~ dbern(PI[i,first[i],j])
      }
    }
    
    # Loop over all other individuals
    for(i in (N[1]+1):N[2]){
      z[i,first[i]] ~ dbern(1)
      S[i,1,first[i]] ~ dunif(xlow[i], xupp[i]) # Prior for the first x coordinate
      S[i,2,first[i]] ~ dunif(ylow[i], yupp[i]) # Prior for the first y coordinate
      # First primary session:
      g[i,first[i],1] <- 0
      for(r in 1:R){ # trap
        D[i,r,first[i]] <- sqrt(pow(S[i,1,first[i]]-X[r,1],2) + pow(S[i,2,first[i]]-X[r,2],2))
        g[i,first[i],r+1] <- exp(-pow(D[i,r,first[i]]/sigma[gr[i]], 2)) # Trap exposure
      }
      G[i,first[i]] <- sum(g[i,first[i],]) # Total trap exposure
      for(j in 1:J[i,first[i]]){
        P[i,j,first[i]] <- 1 - exp(-lambda[gr[i]]*G[i,first[i]]) # Probability of being captured
        PI[i,first[i],j] <- step(H[i,j,first[i]]-2)*(g[i,first[i],H[i,j,first[i]]]/(G[i,first[i]]+ 0.001))*P[i,j,first[i]] + (1.001-step(H[i,j,first[i]]-2))*(1.001-P[i,j,first[i]])
        Ones[i,j,first[i]] ~ dbern(PI[i,first[i],j])
      }
      ## Later primary sessions
      for(k in (first[i]+1):K[i]){ # primary session
        theta[i,k-1] ~ dunif(-3.141593,3.141593) # Prior for dispersal direction
        z[i,k] ~ dbern(Palive[i,k-1])
        Palive[i,k-1] <- z[i,k-1] * phi[gr[i],k-1] # Pr(alive in primary session k)
        d[i,k-1] ~ dexp(dlambda[gr[i]])
        ## assumption that movement btw primary sessions is happening in a random direction at a random distance
        S[i,1,k] <- S[i,1,k-1] + d[i,k-1]*cos(theta[i,k-1])
        S[i,2,k] <- S[i,2,k-1] + d[i,k-1]*sin(theta[i,k-1])
        g[i,k,1] <- 0
        for(r in 1:R){ # trap
          D[i,r,k] <- sqrt(pow(S[i,1,k]-X[r,1],2) + pow(S[i,2,k]-X[r,2],2))  # Squared distance to trap
          g[i,k,r+1] <- exp(-pow(D[i,r,k]/sigma[gr[i]], 2)) # Trap exposure
        }
        G[i,k] <- sum(g[i,k,]) # Total trap exposure
        for(j in 1:J[i,k]){
          P[i,j,k] <- (1 - exp(-lambda[gr[i]]*G[i,k]))*z[i,k] # Probability of being captured
          PI[i,k,j] <- step(H[i,j,k]-2)*(g[i,k,H[i,j,k]]/(G[i,k] + 0.01))*P[i,j,k] + (1.001-step(H[i,j,k]-2))*(1.001-P[i,j,k])
          Ones[i,j,k] ~ dbern(PI[i,k,j])
        } # End j
      } # End k
  
  # Average estimated locations at final survey
      S_[i,1] <- (S[i,1,6])
      S_[i,2] <- (S[i,2,6])
    } # End i
  
    #Derived Parameters
    ## probability of surviving across entire study time period (all 6 primary periods)
    study.r <- phi[1,1]*phi[1,2]*phi[1,3]*phi[1,4]*phi[1,5]
    study.s <- phi[2,1]*phi[2,2]*phi[2,3]*phi[2,4]*phi[2,5]
    
  }
  ", con = "sparcNet_mod.txt")


# Parameters to monitor
params0 = c(
  "sigma",
  "dmean",
  "study.s",
  "study.r",
  "lambda",
  "phi",
  'Phi',
  "S_"
)

#MCMC settings
saved.per.chain <- 500 
nt <- 10 
nb <- 4000 
nc <- 20
na <- 5000 
(ni <- saved.per.chain * nt + na + nb)



modNULL.logit <- jagsUI::jags(
  data0,
  parameters.to.save = params0,
  model.file = "sparcNet_mod.txt",
  parallel = T,
  n.adapt = na,
  n.chains = nc,
  n.thin = nt,
  n.iter = ni,
  n.burnin = nb
)