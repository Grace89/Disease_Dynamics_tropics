# Create a function to simulate the data

data.fn <- function(R = 100, # Number of sites
                    T = 3,   # Number of replicate surveys per site
                    K = 10,   # Number of seasons
                    
                    alpha.lamN = 3,   # Average number of hosts per site the first season
                    alpha.lamI = 3,   # Average number of hosts per site the first season
                    phiN = 0.9,  # Apparent host survival probability of uninfected and infected hosts
                    phiI = 0.7,  # Apparent host survival probability of uninfected and infected hosts
                    gammaN = 3, # Average number of individuals arriving at each site for uninfected and infected hosts
                    gammaI = 2, # Average number of individuals arriving at each site for uninfected and infected hosts
                    pN = 0.8,  # Detection probability during each survey
                    pI = 0.7,  # Detection probability during each survey
                    psi_NI = 0.1,  # Recovery probability
                    psi_IN = 0.1  # Infection probability
){
  
  # Empty matrices to hold the data    
  yN <- yI <- array(NA, dim = c(R, T, K))   # Observed abundance data
  NI <- NN <- array(NA, dim = c(R, K))      # True abundance data
  
  
  #-------- First season 
  NN[,1] <- rpois(n = R, lambda = alpha.lamN)
  NI[,1] <- rpois(n = R, lambda = alpha.lamI)
  
  #------ Second season  
  # Empty matrices to hold latent abundance variables 
  # i.e., number of hosts surviving, arriving, and transitioning
  
  SN <- SI <- GI <- GN <- TrN <- TrI <- array(0, dim = c(R, K-1))  
  
  for(k in 2:K){  
    
    for(i in 1:R){
      
      if(NN[i,k-1]>0){
        SN[i, k-1] <- rbinom(n=1, size=NN[i,k-1], prob=phiN)       
        # Survival of not infecteds
        TrN[i,k-1] <- rbinom(n=1, size=SN[i,k-1], prob=psi_NI)    
        # Getting infected
      }
      if(NI[i,k-1]>0){
        SI[i, k-1] <-  rbinom(n=1, size=NI[i,k-1], prob=phiI)   
        # Survival of infecteds
        TrI[i, k-1] <- rbinom(n=1, size=SI[i,k-1], prob=psi_IN) 
        # Losing infection
      }
      # Recruitment
      GN[i, k-1] <- rpois(1, lambda = gammaN)
      GI[i, k-1] <- rpois(1, lambda = gammaI)
      
    }
    
    # Total
    NI[,k] <-  SI[,k-1] + GI[,k-1] + TrN[,k-1] - TrI[,k-1]
    NN[,k] <-  SN[,k-1] + GN[,k-1] + TrI[,k-1] - TrN[,k-1]
    
  }
  
  # Obervation process
  
  for(i in 1:R){
    
    for(j in 1:T){
      
      for(k in 1:K){
        
        yN[i, j, k] <- rbinom(n = 1, NN[i,k], pN)
        yI[i, j, k] <- rbinom(n = 1, NI[i,k], pI)
        
      }
    }
  }
  
  
  return(list(R = R, T = T, K = K,
              alpha.lamN= alpha.lamN,
              alpha.lamI= alpha.lamI,
              phiN = phiN,
              phiI = phiI,
              gammaN = gammaN,
              gammaI = gammaI,
              psi_NI = psi_NI,
              psi_IN = psi_IN,
              SN = SN,
              SI = SI,
              GN = GN,
              GI = GI,
              TrI = TrI,
              TrN = TrN,
              NN = NN, 
              NI = NI,
              pN = pN, 
              pI = pI, 
              yN = yN,
              yI = yI))
}


sink("model.txt")
cat("
model{

# Priors

#------- NOT Infected
alpha.lamN  ~ dnorm(0,0.01) # Average abundance of uninfected hosts the first season
pN     ~ dunif(0,1) # Detection probability of uninfected hosts
gammaN ~ dnorm(0,0.01)  # Gains rate of uninfected hosts
phiN   ~ dunif(0,1) # Survival probability of uninfected hosts
psi_NI ~ dunif(0,1) # Infection probability (transitioning from uninfected to infected)

#------- Infected
alpha.lamI  ~ dnorm(0,0.01)# Average abundance of infecteds the first season
pI     ~ dunif(0,1) # Detection probability of infected
gammaI ~ dnorm(0,0.01)  # Gains rate of infected
phiI   ~ dunif(0,1) # Survival probability of infected
psi_IN ~ dunif(0,1)  # Recovery probability (transitioning from infected to uninfected)

#------------ Ecological model
#---- First season

for(i in 1:R){
#------ Not infected
  NN[i, 1] ~ dpois(lambdaN[i])
    log(lambdaN[i]) <- alpha.lamN

#----- Infected
  NI[i, 1] ~ dpois(lambdaI[i])  
    log(lambdaI[i]) <- alpha.lamI
}

#------ All other seasons
for(k in 2:K){

  for(i in 1:R){

#------- Nost Infected

    SN[i,k] ~ dbin(phiN, NN[i,k-1])  # Total survivors
    TN[i,k] ~ dbin(psi_NI, SN[i,k])  # Survive, become infected

    GN[i,k] ~ dpois(GaN[i, k])   # Recruits
    log(GaN[i, k]) <- gammaN

#------- Infected

    SI[i,k] ~ dbin(phiI, NI[i,k-1] ) # Infecteds who survive
    TI[i,k] ~ dbin(psi_IN, SI[i,k] ) # Get better, transition to uninfected

    GI[i,k] ~ dpois(GaI[i, k])  # Recruits
      log(GaI[i, k]) <- gammaI

#------- Totals

    NN[i, k]  <- SN[i,k] - TN[i,k] + GN[i,k] + TI[i,k]  
    NI[i, k]  <- SI[i,k] - TI[i,k] + GI[i,k]  + TN[i,k]  

}

}

#------------- Obervation model

for(i in 1:R){

  for(j in 1:T){

    for(k in 1:K){

      yN[i, j, k] ~ dbin(pN, NN[i, k]) 
      yI[i, j, k] ~ dbin(pI, NI[i, k]) 


    }

  }

}


}
", fill = TRUE)
sink()


# Simulate the data
sodata<- data.fn()

# Bundle data
win.data <- list(yN = sodata$yN, 
                 yI = sodata$yI,
                 R = dim(sodata$yN)[1], 
                 T = dim(sodata$yN)[2],
                 K = dim(sodata$yN)[3])

# Initial values
NIst <- apply(sodata$yI, c(1, 3), max, na.rm = TRUE)
NNst <- apply(sodata$yN, c(1, 3), max, na.rm = TRUE)

inits <- function() {list( 
  alpha.lamN = sodata$alpha.lamN,
  pN = sodata$pN,
  phiN = sodata$phiN,
  gammaN = sodata$gammaN,
  psi_NI = sodata$psi_NI,
  
  alpha.lamI = sodata$alpha.lamI,
  pI = sodata$pI,
  phiI = sodata$phiI,
  gammaI = sodata$gammaI,                        
  psi_IN = sodata$psi_IN
  
)}

# Monitor Parameters
params <- c("alpha.lamN", "alpha.lamI",
            "pN", "pI",
            "phiN", "phiI",
            "gammaN", "gammaI",
            "psi_NI", "psi_IN"
)

# MCMC settings
ni <- 2400
nb <- 400
nt <- 10
nc <- 3


library("jagsUI")

out <- jags(win.data, inits, params, "model.txt", 
            n.chains = nc, n.thin = nt, n.iter = ni, 
            n.burnin = nb, parallel = TRUE)


print(out)
plot(out)



# Store the true paramter values
true <- c(log(sodata$alpha.lamN),
          log(sodata$alpha.lamI),
          (sodata$pN),
          (sodata$pI),
          (sodata$phiN),
          (sodata$phiI),
          log(sodata$gammaN),
          log(sodata$gammaI),
          (sodata$psi_NI), 
          (sodata$psi_IN)
)

# Store the parameter names
names <- c("alpha.lamN", "alpha.lamI",
           "pN", "pI",
           "phiN", "phiI",
           "gammaN", "gammaI",
           "psi_NI", "psi_IN"
)

# Extract the model mean
mod.mean <- c(
  out$mean$alpha.lamN,
  out$mean$alpha.lamI,
  out$mean$pN,
  out$mean$pI,
  out$mean$phiN,
  out$mean$phiI,
  out$mean$gammaN,
  out$mean$gammaI,
  out$mean$psi_NI,
  out$mean$psi_IN)

# Extract the lower credible interval
mod.q2.5 <- c(
  out$q2.5$alpha.lamN,
  out$q2.5$alpha.lamI,
  out$q2.5$pN,
  out$q2.5$pI,
  out$q2.5$phiN,
  out$q2.5$phiI,
  out$q2.5$gammaN,
  out$q2.5$gammaI,
  out$q2.5$psi_NI,
  out$q2.5$psi_IN)

# Extract the upper credible interval
mod.q97.5 <- c(
  out$q97.5$alpha.lamN,
  out$q97.5$alpha.lamI,
  out$q97.5$pN,
  out$q97.5$pI,
  out$q97.5$phiN,
  out$q97.5$phiI,
  out$q97.5$gammaN,
  out$q97.5$gammaI,
  out$q97.5$psi_NI,
  out$q97.5$psi_IN)

# Combine names, truth, mean, lower CI, and upper CI
dat <- data.frame(names = names, 
                  true = true, 
                  mod.mean = mod.mean, 
                  mod.q2.5 = mod.q2.5, 
                  mod.q97.5 = mod.q97.5)

# Load library
library(ggplot2)

# Indicate colors for the plot
cols <- c("Truth" = "red", "Estimated" = "black")

# Plot
ggplot(dat, aes(x= names, y=mod.mean))+ 
  geom_linerange(size = 1, aes(ymin=mod.q2.5, ymax=mod.q97.5)) +
  geom_point(size = 3, aes(x = names, y = mod.mean, col = "Estimated")) +
  geom_point(size = 3, aes(x = names, y = true, col = "Truth")) +
  scale_colour_manual("Values", values=cols)+
  geom_hline(yintercept = 0, lty=2) +
  coord_flip() + ylab('Parameter estimates') +
  xlab("Parameter names") +
  theme_bw()+ 
  theme(axis.text.x = element_text(size = 17, color = "black"), 
        axis.text.y = element_text(size = 17, color = "black"), 
        axis.title.y = element_text(size = 17, color = "black"), 
        axis.title.x =element_text(size = 17, color = "black"),
        legend.title =element_text(size = 17, color = "black"),
        legend.text =element_text(size = 17, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 


