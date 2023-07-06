 rm(list=ls()) #do you want to nuke the environ?
library(R2jags)
library(coda)
 library(devEMF)


### BRING IN ANCHOVY SSB DATA #################################################

# First we are going to use SSB estimates from Thayer et al 2017. Note that the
# annual CVs for the SSB estimates are not actually CVs. We know this because the CVs they
# provide from the McCall paper don't match the error bars shown in McCall. Through trial and error
# it appears the values in the CV column are actaully SE of the log SSB. That is, the values directly
# match the SE error bars provided in McCall. 

# We have digitized these data (log SSB and SE), and bring them in here:


# As a lowbrow means of accounting for the SE in these data (without specific observations)
# we generate 100 random log SSB data points per year that have the exact annual mean and SD 
# provided in the Thayer et al 2017 table. 

SSB<-read.csv("AnchovySSB_Thayer2017_short.csv", header = T)

# Now I'm going to change the years so that the earliest year equals 1...
max.yr.ab<-max(SSB$Year); min.yr.ab<-min(SSB$Year)
SSB$Year<-SSB$Year-min(SSB$Year) +1

# Get the list of years data are assigned to
years.ab<-SSB$Year

#get the number of years the data set spans
N.yr.ab<-max(years.ab)

# Now let's make some "data" #########################
sampsize<-20
weights.ab<-matrix(data=rnorm(sampsize*length(years.ab)),nrow=sampsize,ncol=length(years.ab)) #data matrix
for (i in 1:length(years.ab)){
  weights.ab[,i] <-weights.ab[,i]*SSB$SD_ssb[i] + SSB$log_ssb[i]
}
boxplot(exp(weights.ab),outline=FALSE)
w<-cbind(as.vector(weights.ab),rep(years.ab, each=sampsize)) #stack the data, add year

# Get the length of the data set
N.ver.ab <-dim(w)[1]


# Specify model in BUGS language ##########################
sink("ssm.jags.ab")
cat("
model { 
# Priors and constraints
N.est[1] ~ dunif(0, 4000)            # Prior for initial yearly max size
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
   weights[i,1] ~ dnorm(N.est[weights[i,2]], tau.obs)
   }
}",fill = TRUE)
sink()

# Bundle data
jags.data <- list(weights=w,N.ver=N.ver.ab, N.yr=N.yr.ab)


# Parameters monitored
parameters <- c("sigma2.obs", "sigma2.proc", "N.est")

# MCMC settings
ni <- 25000
nt <- 3
nb <- 10000
nc <- 3


# Call JAGS from R ################
ssm2 <- jags(jags.data, inits=NULL, parameters, "ssm.jags.ab", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
attach.jags(ssm2)
fitted_SSB<-apply(N.est,2,mean)
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
  plot(0,0,ylim=c(m1,m2),xlim=c(min.yr,max.yr),ylab="SSB (1000 metric tons)",xlab="year",las=1,col="black",type="l",lwd=2,frame=FALSE,axes=FALSE)
  axis(2,las=1)
  axis(1,at=seq(min.yr,max.yr,5), labels=seq(min.yr,max.yr,5))
  #axis(1,at=0:N.yr,labels=rep("",N.yr+1),tcl=-0.25)
  polygon(x=c(min.yr:max.yr,max.yr:min.yr),y=c(lower,upper[N.yr:1]),col="gray90",border="gray90")
  #points(years+(min.yr-1),weights,pch = ".")
  points(seq(min.yr,max.yr),fitted,type="l",col="red",lwd=2)
  
}

graph.ssm.SSB(exp(N.est),N.yr.ab,w,years.ab,max.yr.ab,min.yr.ab)

emf("SSB_StateSpace.emf",width = 14, height = 7,emfPlus = FALSE)
graph.ssm.SSB(exp(N.est),N.yr.ab,w,years.ab,max.yr.ab,min.yr.ab)
dev.off()

# Now write results to file
fitted<-lower95<-upper95<-numeric()
for (i in 1:N.yr.ab){
  fitted[i]<-mean(N.est[,i])
  lower95[i]<-quantile(N.est[,i],0.025)
  upper95[i]<-quantile(N.est[,i],0.975)
}
Year<-seq(from = min.yr.ab, to = max.yr.ab, by = 1) 
SSB_SS_MAP<- fitted
SSB_SS_lower95<- lower95
SSB_SS_upper95<- upper95

SSB.SS.estimates<-cbind(Year,SSB_SS_MAP,SSB_SS_lower95,SSB_SS_upper95)
write.csv(SSB.SS.estimates, file = "SSB.SS.estimates.csv")

# lags in X (FCL) that predict Y (SSB):
ccf(fitted_FCL[1:46],fitted_SSB) # dang... that's rad.
ccfvalues<-ccf(fitted_FCL[1:46],fitted_SSB)
# Thus, FCL from 2 years ago best predicts SSB in current year (and it predicts it pretty darn well)
plot(fitted_FCL[1:46],exp(fitted_SSB))



