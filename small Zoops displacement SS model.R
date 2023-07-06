rm(list=ls()) #do you want to nuke the environ?
library(R2jags)
library(coda)
library(devEMF)


### BRING IN DATA ####################################################################

# let's bring the small fraction zooplankton biomass  
zoops <- read.csv("Small zoopl. displ. vol. 20210921.csv", header = T)
str(zoops) #peep the deets of the dat

# 5 of 5345 observations are 0 (which is not really possible). I remove these before logging for convenience
zoops<-zoops[-(which(zoops$Sml_PVolC3==0)),]

# Now I'm going to change the years so that the earliest year equals 1...
max.yr<-max(zoops$Year); min.yr<-min(zoops$Year)
Year<-seq(from = min.yr, to = max.yr) 
zoops$Year<-zoops$Year-min(zoops$Year) +1

# Get the list of years data are assigned to
years<-as.matrix(zoops$Year)

#get the number of years the data set spans
N.yr<-max(years)

# Get the length of the data set
N.ver <-dim(zoops)[1]

#get small fraction displacement volume (log)
w<-as.matrix(log(zoops$Sml_PVolC3))


# Specify model in BUGS language ##########################
sink("ssm.jags.z")
cat("
model { 
# Priors and constraints
N.est[1] ~ dunif(.1, 9)            # Prior for initial year zoops
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
   weights[i,1] ~ dnorm(N.est[years[i,1]], tau.obs)
   }
}",fill = TRUE)
sink()

# Bundle data
jags.data <- list(weights=w,N.ver=N.ver, N.yr=N.yr, years=years)


# Parameters monitored
parameters <- c("sigma2.obs", "sigma2.proc", "N.est")

# MCMC settings
ni <- 25000
nt <- 3
nb <- 10000
nc <- 3


# Call JAGS from R ################
ssm2 <- jags(jags.data, inits=NULL, parameters, "ssm.jags.z", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
attach.jags(ssm2)
fitted_zoop<-apply(N.est,2,mean)
# Graph results ###############
graph.ssm.SSB<-function(N.est,N.yr,w,years,max.yr, min.yr){
  fitted<-lower<-upper<-numeric()
  
  for (i in 1:N.yr){
    fitted[i]<-mean(N.est[,i])
    lower[i]<-quantile(N.est[,i],0.1)
    upper[i]<-quantile(N.est[,i],0.9)
  }
  m1<-min(c(fitted,lower))
  m2<-max(c(fitted,upper)) 
  par(mar=c(4.5,4,1,1),cex=1.2)
  plot(0,0,ylim=c(m1,m2),xlim=c(min.yr,max.yr),ylab="Sml_PVol",xlab="year",las=1,col="black",type="l",lwd=2,frame=FALSE,axes=FALSE)
  axis(2,las=1)
  axis(1,at=seq(min.yr,max.yr,5), labels=seq(min.yr,max.yr,5))
  #axis(1,at=0:N.yr,labels=rep("",N.yr+1),tcl=-0.25)
  polygon(x=c(min.yr:max.yr,max.yr:min.yr),y=c(lower,upper[N.yr:1]),col="gray90",border="gray90")
  #points(years+(min.yr-1),weights,pch = ".")
  points(seq(min.yr,max.yr),fitted,type="l",col="red",lwd=2)
  
}
graph.ssm.SSB(exp(N.est),N.yr,w,years,max.yr,min.yr)
emf("zoop_Displacement_StateSpace.emf",width = 14, height = 7,emfPlus = FALSE)
graph.ssm.SSB(exp(N.est),N.yr,w,years,max.yr,min.yr)
dev.off()


# Now write results to file
fitted<-lower95<-upper95<-numeric()
for (i in 1:N.yr){
  fitted[i]<-mean(exp(N.est[,i]))
  lower95[i]<-quantile(exp(N.est[,i]),0.025)
  upper95[i]<-quantile(exp(N.est[,i]),0.975)
}

zoop_Displacement_SS_MAP<- fitted
zoop_Displacement_SS_lower95<- lower95
zoop_Displacement_SS_upper95<- upper95

zoop_Displacement.SS.estimates<-cbind(Year,zoop_Displacement_SS_MAP,zoop_Displacement_SS_lower95,zoop_Displacement_SS_upper95)
write.csv(zoop_Displacement.SS.estimates, file = "zoop_Displacement.SS.estimates.csv")



