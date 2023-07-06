rm(list=ls()) #do you want to nuke the environ?
library(R2jags)
library(coda)
library(devEMF)

### BRING IN COPEPOD DATA ####################################################################
# First let's read in the raw data file (observation records in flat format)
R <- read.csv("copepod time series data.csv", header = T)
str(R) #peep the deets of the dat

# Remove first two years to match other analyses as part of the project
R <- subset(R,year > 1959 )

# Get the length of the data set
N.ver <-dim(R)[1]

# Get the list of years data are assigned to
years<-R$year-min(R$year) +1

#get the number of years the data set spans
N.yr<-max(years)

# Get the log copepod data
weights<-log(R$Metridia.cf..pacifica)

# Specify model in BUGS language
sink("ssm.jags")
cat("
model { 
# Priors and constraints
N.est[1] ~ dunif(0, 20)            # Prior for initial yearly max size
sigma.proc ~ dunif(0, 5)           # Prior for sd of state process
sigma2.proc <- pow(sigma.proc, 2)
tau.proc <- pow(sigma.proc, -2)
sigma.obs ~ dunif(0, 5)           # Prior for sd of observation process
sigma2.obs <- pow(sigma.obs, 2)
tau.obs <- pow(sigma.obs, -2)

# Likelihood
# State process
for (t in 1:(N.yr-1)){
   mu[t] ~ dnorm(0, tau.proc) 
   N.est[t+1] <- N.est[t] + mu[t] 
   }

# Observation process
for (i in 1:N.ver) {
   weights[i] ~ dnorm(N.est[years[i]], tau.obs)
   }
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(weights=weights,years=years,N.ver=N.ver, N.yr=N.yr)


# Parameters monitored
parameters <- c("sigma2.obs", "sigma2.proc", "N.est")

# MCMC settings
ni <- 25000
nt <- 3
nb <- 10000
nc <- 3


# Specify model in BUGS language
library(R2jags)
library(coda)

# Call JAGS from R (BRT <1 min)
ssm <- jags(jags.data, inits=NULL, parameters, "ssm.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

attach.jags(ssm)
graph.ssm<-function(N.est,N.yr,weights,years,max.yr, min.yr){
  fitted<-lower<-upper<-numeric()
  
  for (i in 1:N.yr){
    fitted[i]<-mean(N.est[,i])
    lower[i]<-quantile(N.est[,i],0.025)
    upper[i]<-quantile(N.est[,i],0.975)
  }
  #m1<-min(c(weights,fitted,lower))
  #m2<-max(c(weights,fitted,upper))
  m1<-min(c(fitted,lower))
  m2<-max(c(fitted,upper)) 
  par(mar=c(4.5,4,1,1),cex=1.2)
  plot(0,0,ylim=c(m1,m2),xlim=c(0.5,N.yr),ylab="log(Metridia)",xlab="year",las=1,col="black",type="l",lwd=2,frame=FALSE,axes=FALSE)
  axis(2,las=1)
  axis(1,at=seq(0,N.yr,5), labels=seq(0,N.yr,5))
  axis(1,at=0:N.yr,labels=rep("",N.yr+1),tcl=-0.25)
  polygon(x=c(1:N.yr,N.yr:1),y=c(lower,upper[N.yr:1]),col="gray90",border="gray90")
  #points(years,weights,pch = ".")
  points(fitted,type="l",col="red",lwd=2)
  
}
graph.ssm(N.est,N.yr,weights,years,Y.max) #output to screen
#Print figure to file
library(devEMF)
emf("Metridia.SS.graph.emf",width = 14, height = 7,emfPlus = FALSE)
graph.ssm(N.est,N.yr,weights,years,max.yr,min.yr)
dev.off()

# Now write results to file
fitted<-lower95<-upper95<-numeric()
for (i in 1:N.yr){
  fitted[i]<-mean(N.est[,i])
  lower95[i]<-quantile(N.est[,i],0.025)
  upper95[i]<-quantile(N.est[,i],0.975)
}
Year<-seq(from = min(R$year), to = max(R$year), by = 1) # extent of copepod data
COP_SS_MAP<- fitted
COP_SS_lower95<- lower95
COP_SS_upper95<- upper95


Metridia.SS.estimates<-cbind(Year,COP_SS_MAP,COP_SS_lower95,COP_SS_upper95)
write.csv(Metridia.SS.estimates, file = "Metridia.SS.estimates.csv")

