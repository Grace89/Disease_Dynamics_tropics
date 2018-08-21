
R = 100 # Number of sites
T = 3   # Number of replicate surveys per site
K = 10   # Number of seasons
S = 40   # Number of species

alpha.lamN = 5   # Average number of hosts per site the first season (log scale)
alpha.lamI = 4   # Average number of hosts per site the first season (log scale)
phiN = 1  # Apparent host survival probability of uninfected hosts (logit scale)
phiI = 0.5  # Apparent host survival probability of infected hosts (logit scale)
gammaN = 5 # Average number of individuals arriving at each site for uninfected and infected hosts (log scale)
gammaI = 4 # Average number of individuals arriving at each site for uninfected and infected hosts (log scale)
pN = 2  # Detection probability during each survey (logit scale)
pI = 1  # Detection probability during each survey (logit scale)
psi_NI = -2 # Recovery probability (logit scale)
psi_IN = -2  # Infection probability (logit scale)
# Create a function to simulate the data

data.fn <- function(R = 100, # Number of sites
                    T = 3,   # Number of replicate surveys per site
                    K = 10,   # Number of seasons
                    S = 40,   # Number of species
                    
                    alpha.lamN = 5,   # Average number of hosts per site the first season (log scale)
                    alpha.lamI = 4,   # Average number of hosts per site the first season (log scale)
                    phiN = 1,  # Apparent host survival probability of uninfected hosts (logit scale)
                    phiI = 0.5,  # Apparent host survival probability of infected hosts (logit scale)
                    gammaN = 5, # Average number of individuals arriving at each site for uninfected and infected hosts (log scale)
                    gammaI = 4, # Average number of individuals arriving at each site for uninfected and infected hosts (log scale)
                    pN = 2,  # Detection probability during each survey (logit scale)
                    pI = 1,  # Detection probability during each survey (logit scale)
                    psi_NI = -2,  # Recovery probability (logit scale)
                    psi_IN = -2  # Infection probability (logit scale)
){

# Empty matrices to hold the data    
yN <- yI <- array(NA, dim = c(R, T, K, S))   # Observed abundance data
NI <- NN <- array(NA, dim = c(R, K, S))      # True abundance data

# Create histograms of average abundance, survival, detection, transmission, and recover for each species in each disease state

# Average abundance the first season
LambdaN <- rnorm(S, mean = alpha.lamN, sd = 1)
LambdaI <- rnorm(S, mean = alpha.lamI, sd = 1)

# Survival
PhiN <- rnorm(S, mean = phiN, sd = 1)
PhiI <- rnorm(S, mean = phiI, sd = 1)

# Recruitment
GammaN <- rnorm(S, mean = gammaN, sd = 1)
GammaI <- rnorm(S, mean = gammaI, sd = 1)

# Detection
PN <- rnorm(S, mean = pN, sd = 1)
PI <- rnorm(S, mean = pI, sd = 1)

# Transitions
Psi_NI <- rnorm(S, mean = psi_NI, sd = 1)
Psi_IN <- rnorm(S, mean = psi_IN, sd = 1)

#--- Plot
par(mfrow = c(5, 2), mar = c(4,4,4,4))
hist(LambdaN, main = "")
abline(v = alpha.lamN, col = "red", lwd = 2)
hist(LambdaI, main = "")
abline(v = alpha.lamI, col = "red", lwd = 2)

hist(plogis(PhiN), main = "")
abline(v = plogis(phiN), col = "red", lwd = 2)
hist(plogis(PhiI), main = "")
abline(v = plogis(phiI), col = "red", lwd = 2)

hist(GammaN, main = "")
abline(v = gammaN, col = "red", lwd = 2)
hist(GammaI, main = "")
abline(v = gammaI, col = "red", lwd = 2)

hist(plogis(PN), main = "")
abline(v = plogis(pN), col = "red", lwd = 2)
hist(plogis(PI), main = "")
abline(v = plogis(pI), col = "red", lwd = 2)

hist(plogis(Psi_NI), main = "")
abline(v = plogis(psi_NI), col = "red", lwd = 2)
hist(plogis(Psi_IN), main = "")
abline(v = plogis(psi_IN), col = "red", lwd = 2)

#-------- First season 
for(s in 1:S){
  NN[,1,s] <- rpois(n = R, lambda = LambdaN[s])
  NI[,1,s] <- rpois(n = R, lambda = LambdaI[s])
}
#------ Second season  
# Empty matrices to hold latent abundance variables 
# i.e., number of hosts surviving, arriving, and transitioning

SN <- SI <- GI <- GN <- TrN <- TrI <- array(0, dim = c(R, K-1, S))  

for(k in 2:K){  
  
  for(i in 1:R){
    
    for(s in 1:S){
      
    if(NN[i,k-1,s]>0){
      SN[i, k-1,s] <- rbinom(n=1, size=NN[i,k-1,s], prob=plogis(PhiN[s]))       
      # Survival of not infecteds
      TrN[i,k-1,s] <- rbinom(n=1, size=SN[i,k-1,s], prob=plogis(Psi_NI[s]))    
      # Getting infected
    }
    if(NI[i,k-1,s]>0){
      SI[i, k-1,s] <-  rbinom(n=1, size=NI[i,k-1,s], prob=plogis(PhiI[s]))   
      # Survival of infecteds
      TrI[i, k-1,s] <- rbinom(n=1, size=SI[i,k-1,s], prob=plogis(Psi_IN[s]))
      # Losing infection
    }
    # Recruitment
    GN[i, k-1,s] <- rpois(1, lambda = GammaN[s])
    GI[i, k-1,s] <- rpois(1, lambda = GammaI[s])

    # Total
    NI[,k,s] <-  SI[,k-1,s] + GI[,k-1,s] + TrN[,k-1,s] - TrI[,k-1,s]
    NN[,k,s] <-  SN[,k-1,s] + GN[,k-1,s] + TrI[,k-1,s] - TrN[,k-1,s]
    
    }
  
  }

}

# Obervation process

for(i in 1:R){
  
  for(j in 1:T){
    
    for(k in 1:K){
      
      for(s in 1:S){
        
      yN[i, j, k, s] <- rbinom(n = 1, NN[i,k,s], plogis(PN[s]))
      yI[i, j, k, s] <- rbinom(n = 1, NI[i,k,s], plogis(PI[s]))
      
      }
    }
  }
}

return(list(R = R, T = T, K = K,
            LambdaN = LambdaN,
            LambdaI = LambdaI,
            PhiN = PhiN,
            PhiI = PhiI,
            GammaN = GammaN,
            GammaI = GammaI,
            PN = PN,
            PI = PI,
            Psi_NI = Psi_NI,
            Psi_IN = Psi_IN,
            alpha.lamN= alpha.lamN,
            alpha.lamI= alpha.lamI,
            phiN = phiN,
            phiI = phiI,
            PhiN = PhiN,
            PhiI = PhiI,
            gammaN = gammaN,
            gammaI = gammaI,
            GammaN = GammaN,
            GammaI = GammaI,
            psi_NI = psi_NI,
            psi_IN = psi_IN,
            Psi_NI = Psi_NI,
            Psi_IN = Psi_IN,
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
            PN = PN, 
            PI = PI,
            yN = yN,
            yI = yI))
}

# Simulate the data
sodata <- data.fn()


sink("model.txt")
cat("
model{

# Priors

#------- NOT Infected

for(s in 1:S){
  alpha.lamN[s]  ~ dnorm(alphaN, alphaN.tau) 
  gammaN[s]  ~ dnorm(gamN, gamN.tau) 
  pN[s] ~ dnorm(mu.pN, tau.pN)
  phiN[s] ~ dnorm(mu.phiN, tau.phiN)
  psi_NI[s] ~ dnorm(mu.psi_NI, tau.psi_NI)
}
    

# Log scale
alphaN ~ dnorm(0, 0.01)
alphaN.tau <- 1/(alphaN.sd * alphaN.sd)
alphaN.sd ~ dgamma(0.01, 0.01)

gamN ~ dnorm(0, 0.01)
gamN.tau <- 1/(gamN.sd * gamN.sd)
gamN.sd ~ dgamma(0.01, 0.01)

# Logit scale
pN.mean  ~ dunif(0,1) 
mu.pN <- log(pN.mean) - log(1-pN.mean)
tau.pN <- 1/(sd.pN * sd.pN)
sd.pN ~ dgamma(0.01, 0.01)

phiN.mean  ~ dunif(0,1) 
mu.phiN <- log(phiN.mean) - log(1-phiN.mean)
tau.phiN <- 1/(sd.phiN * sd.phiN)
sd.phiN ~ dgamma(0.01, 0.01)

psi_NI.mean  ~ dunif(0,1) 
mu.psi_NI <- log(psi_NI.mean) - log(1-psi_NI.mean)
tau.psi_NI <- 1/(sd.psi_NI * sd.psi_NI)
sd.psi_NI ~ dgamma(0.01, 0.01)

#------- Infected

for(s in 1:S){
  alpha.lamI[s]  ~ dnorm(alphaI, alphaI.tau) 
  gammaI[s]  ~ dnorm(gamI, gamI.tau) 
  pI[s] ~ dnorm(mu.pI, tau.pI)
  phiI[s] ~ dnorm(mu.phiI, tau.phiI)
  psi_IN[s] ~ dnorm(mu.psi_IN, tau.psi_IN)
}


# Log scale
alphaI ~ dnorm(0, 0.01)
alphaI.tau <- 1/(alphaI.sd * alphaI.sd)
alphaI.sd ~ dgamma(0.01, 0.01)

gamI ~ dnorm(0, 0.01)
gamI.tau <- 1/(gamI.sd * gamI.sd)
gamI.sd ~ dgamma(0.01, 0.01)

# Logit scale
pI.mean  ~ dunif(0,1) 
mu.pI <- log(pI.mean) - log(1-pI.mean)
tau.pI <- 1/(sd.pI * sd.pI)
sd.pI ~ dgamma(0.01, 0.01)

phiI.mean  ~ dunif(0,1) 
mu.phiI <- log(phiI.mean) - log(1-phiI.mean)
tau.phiI <- 1/(sd.phiI * sd.phiI)
sd.phiI ~ dgamma(0.01, 0.01)

psi_IN.mean  ~ dunif(0,1) 
mu.psi_IN <- log(psi_IN.mean) - log(1-psi_IN.mean)
tau.psi_IN <- 1/(sd.psi_IN * sd.psi_IN)
sd.psi_IN ~ dgamma(0.01, 0.01)

#------------ Ecological model
#---- First season

for(s in 1:S){

  for(i in 1:R){
#------ Not infected
  NN[i, 1, s] ~ dpois(lambdaN[i,s])
    log(lambdaN[i,s]) <- alpha.lamN[s]

#----- Infected
  NI[i, 1, s] ~ dpois(lambdaI[i,s])  
    log(lambdaI[i,s]) <- alpha.lamI[s]
  }

#------ All other seasons
for(k in 2:K){

  for(i in 1:R){

#------- Nost Infected

    SN[i,k-1,s] ~ dbin(phiN[s], NN[i,k-1,s])  # Total survivors
    TN[i,k-1,s] ~ dbin(psi_NI[s], SN[i,k-1,s])  # Survive, become infected

    GN[i,k-1,s] ~ dpois(GaN[i, k-1,s])   # Recruits
    log(GaN[i, k-1, s]) <- gammaN[s]

#------- Infected

    SI[i,k-1,s] ~ dbin(phiI[s], NI[i,k-1,s] ) # Infecteds who survive
    TI[i,k-1,s] ~ dbin(psi_IN[s], SI[i,k-1,s] ) # Get better, transition to uninfected

    GI[i,k-1,s] ~ dpois(GaI[i, k-1, s])  # Recruits
      log(GaI[i, k-1, s]) <- gammaI[s]

#------- Totals

    NN[i, k, s]  <- SN[i,k-1,s] - TN[i,k-1,s] + GN[i,k-1,s] + TI[i,k-1,s]  
    NI[i, k, s]  <- SI[i,k-1,s] - TI[i,k-1,s] + GI[i,k-1,s] + TN[i,k-1,s]  

    }

  }

}

#------------- Obervation model

for(s in 1:S){
  for(i in 1:R){
   for(j in 1:T){
    for(k in 1:K){
        yN[i, j, k, s] ~ dbin(pN[s], NN[i, k, s]) 
        yI[i, j, k, s] ~ dbin(pI[s], NI[i, k, s]) 
      }
    }
  }
}


}
", fill = TRUE)
sink()

# Bundle data
win.data <- list(yN = sodata$yN, 
                 yI = sodata$yI,
                 R = dim(sodata$yN)[1], 
                 T = dim(sodata$yN)[2],
                 K = dim(sodata$yN)[3],
                 S = dim(sodata$yN)[4])

# Initial values
NNst <- sodata$NN
NNst[, -1, ] <- NA

NIst <- sodata$NI
NIst[, -1, ] <- NA

inits <- function() {list( 
  alphaN = log(sodata$alpha.lamN),
  pN.mean = plogis(sodata$pN),
  phiN.mean = plogis(sodata$phiN),
  gamN = log(sodata$gammaN),
  psi_NI.mean = plogis(sodata$psi_NI),
  
  alphaI = log(sodata$alpha.lamI),
  pI.mean = plogis(sodata$pI),
  phiI.mean = plogis(sodata$phiI),
  gamI = log(sodata$gammaI),                        
  psi_IN.mean = plogis(sodata$psi_IN),
  
  alpha.lamN = sodata$LambdaN,
  gammaN = sodata$GammaN,
  pN = sodata$PN,
  phiN = sodata$PhiN,
  psi_NI = sodata$Psi_NI,
  
  alpha.lamI = sodata$LambdaI,
  gammaI = sodata$GammaI,
  pI = sodata$PI,
  phiI = sodata$PhiI,
  psi_IN = sodata$Psi_IN,
  
  SN = sodata$SN,
  SI = sodata$SI,
  GN = sodata$GN,
  GI = sodata$GI, 
  TN = sodata$TrN,
  TI = sodata$TrI,
  NN = NNst,
  NI = NIst,
  
  alphaN.sd = 1,
  gamN.sd = 1,
  sd.pN = 1,
  sd.phiN = 1,
  sd.psi_NI = 1,
  alphaI.sd = 1,
  gamI.sd = 1,
  sd.pI = 1,
  sd.phiI = 1,
  sd.psi_IN = 1
  
)}

# Monitor Parameters
params <- c("alphaN", "alphaN.sd",
            "gamN", "gamN.sd",
            "pN.mean", "sd.pN",
            "phiN.mean", "sd.phiN",
            "psi_NI.mean", "sd.psi_NI",
            "alphaI", "alphaI.sd",
            "gamI", "gamI.sd",
            "pI.mean", "sd.pI",
            "phiI.mean", "sd.phiI",
            "psi_IN.mean", "sd.psi_IN"
)

# MCMC settings
ni <- 10
na <- 1
nb <- 1
nt <- 1
nc <- 3


library("jagsUI")

out <- jags(win.data, inits, params, "model.txt", 
            n.chains = nc, n.thin = nt, n.iter = ni, 
            n.burnin = nb, n.adapt = na, parallel = TRUE)

sodata$yN[,,,1]


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


