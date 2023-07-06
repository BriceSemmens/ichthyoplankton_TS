 rm(list=ls()) #do you want to nuke the environ?
library(R2jags)
library(coda)
 library(devEMF)

 
### BRING IN FCL DATA ####################################################################
# First let's read in the raw data file (observation records in flat format)
R.full <- read.csv("Raw FCL data 20211208.csv", header = T)
str(R.full) #peep the deets of the dat

# Pull only anchovy 
R<- R.full[R.full[, "Species"] == 'Anchovy',]

# let's only use the larger ones (explore the impact of shift to smaller animals in later years)
#R<- R.full[R.full[, "SL_mm"] > 17.9,]
# let's only use the larger ones (explore the impact of shift to smaller animals in later years)
#R<- R[R[, "SL_mm"] < 21.1,]
boxplot(R$TP_Trp_Scr~R$Year,outline=FALSE)

# Now I'm going to change the years so that the earliest year equals 1...
max.yr<-max(R$Year); min.yr<-min(R$Year)
years<-R$Year-min(R$Year) +1
unique(R$Year) #note that at least one year in the series is missing

# Get the length of the data set
N.ver <-dim(R)[1]

#get the number of years the data set spans
N.yr<-max(years)

# Get the TL data
#weights<-R$TP_Trp_Scr
#weights<-R$Energy_20
#weights<-R$TP_Glu_Phe
weights<-R$Phe_d15N 

# Specify model in BUGS language #####################

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
}",fill = TRUE)
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

# Call JAGS from R ###################
ssm <- jags(jags.data, inits=NULL, parameters, "ssm.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
attach.jags(ssm)
fitted_FCL<-apply(N.est,2,mean)

# Graph results #####################

graph.ssm.FCL<-function(N.est,N.yr,weights,years,max.yr, min.yr){
  fitted<-lower<-upper<-numeric()
 
  for (i in 1:N.yr){
    fitted[i]<-mean(N.est[,i])
    lower[i]<-quantile(N.est[,i],0.27)
    upper[i]<-quantile(N.est[,i],0.75)
  }
  m1<-min(c(fitted,lower))
  m2<-max(c(fitted,upper)) 
  par(mar=c(4.5,4,1,1),cex=1.2)
  plot(0,0,ylim=c(m1,m2),xlim=c(min.yr,max.yr),ylab="Phe_d15N",xlab="year",las=1,col="black",type="l",lwd=2,frame=FALSE,axes=FALSE)
  axis(2,las=1)
  axis(1,at=seq(min.yr,max.yr,5), labels=seq(min.yr,max.yr,5))
  #axis(1,at=0:N.yr,labels=rep("",N.yr+1),tcl=-0.25)
  polygon(x=c(min.yr:max.yr,max.yr:min.yr),y=c(lower,upper[N.yr:1]),col="gray90",border="gray90")
  #points(years+(min.yr-1),weights,pch = ".")
  points(seq(min.yr,max.yr),fitted,type="l",col="red",lwd=2)

}
graph.ssm.FCL(N.est,N.yr,weights,years,max.yr,min.yr)

emf("Phe_d15N_StateSpace.emf",width = 14, height = 7,emfPlus = FALSE)
graph.ssm.FCL(N.est,N.yr,weights,years,max.yr,min.yr)
dev.off()

# Now write results to file
fitted<-lower95<-upper95<-numeric()
for (i in 1:N.yr){
  fitted[i]<-mean(N.est[,i])
  lower95[i]<-quantile(N.est[,i],0.025)
  upper95[i]<-quantile(N.est[,i],0.975)
}
Year<-seq(from = min(R$Year), to = max(R$Year), by = 1) 
Phe_d15N_SS_MAP<- fitted
Phe_d15N_SS_lower95<- lower95
Phe_d15N_SS_upper95<- upper95

Phe_d15N.SS.estimates<-cbind(Year,Phe_d15N_SS_MAP,Phe_d15N_SS_lower95,Phe_d15N_SS_upper95)
write.csv(Phe_d15N.SS.estimates, file = "Phe_d15N.SS.estimates.csv")

