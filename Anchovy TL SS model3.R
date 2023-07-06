 rm(list=ls()) #do you want to nuke the environ?
library(R2jags)
library(coda)
 library(devEMF)

 # Make a data frame to store drived time series data
 Year<-seq(from = 1951, to = 2009, by = 1) # extent of CalCOFI larval data
 FCL_raw_means<- rep(NaN,length(Year))
 FCL_SS_means<- rep(NaN,length(Year))
 SSB_raw_means<- rep(NaN,length(Year))
 SSB_SS_means<- rep(NaN,length(Year))
 larvae_10mm_plus<- rep(NaN,length(Year))
 larvae_10mm_minus<- rep(NaN,length(Year))
 larvae_all<- rep(NaN,length(Year))
 DF<-data.frame(Year,FCL_raw_means,FCL_SS_means,SSB_raw_means,SSB_SS_means,larvae_10mm_plus,larvae_10mm_minus)
 
 
### BRING IN FCL DATA ####################################################################
# First let's read in the raw data file (observation records in flat format)
R.full <- read.csv("Raw FCL data 20211208.csv", header = T)
str(R.full) #peep the deets of the dat

# Pull only anchovy 
R<- R.full[R.full[, "Species"] == 'Anchovy',]
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
weights<-R$Energy_20
#weights<-R$TP_Glu_Phe
#weights<-R$Phe_d15N 

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
  plot(0,0,ylim=c(m1,m2),xlim=c(min.yr,max.yr),ylab="Transfer Efficiency (%)",xlab="year",las=1,col="black",type="l",lwd=2,frame=FALSE,axes=FALSE)
  axis(2,las=1)
  axis(1,at=seq(min.yr,max.yr,5), labels=seq(min.yr,max.yr,5))
  #axis(1,at=0:N.yr,labels=rep("",N.yr+1),tcl=-0.25)
  polygon(x=c(min.yr:max.yr,max.yr:min.yr),y=c(lower,upper[N.yr:1]),col="gray90",border="gray90")
  #points(years+(min.yr-1),weights,pch = ".")
  points(seq(min.yr,max.yr),fitted,type="l",col="red",lwd=2)

}
emf("Energy_StateSpace.emf",width = 14, height = 7,emfPlus = FALSE)
graph.ssm.FCL(N.est,N.yr,weights,years,max.yr,min.yr)
dev.off()

#print posteriors
trans.Eff.Mean.fitted<-apply(N.est,2,mean)
trans.Eff.SD.fitted<-apply(N.est,2,sd)
year.span<-seq(min.yr,max.yr,1)
Trans.Eff.Means<-cbind(trans.Eff.Mean.fitted,trans.Eff.SD.fitted,year.span)
write.csv(Trans.Eff.Means, file = "Trans.Eff.Means.csv")

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
N.est[1] ~ dunif(0, 4)            # Prior for initial yearly max size
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
  plot(0,0,ylim=c(m1,m2),xlim=c(min.yr,max.yr),ylab="SSB",xlab="year",las=1,col="black",type="l",lwd=2,frame=FALSE,axes=FALSE)
  axis(2,las=1)
  axis(1,at=seq(min.yr,max.yr,5), labels=seq(min.yr,max.yr,5))
  #axis(1,at=0:N.yr,labels=rep("",N.yr+1),tcl=-0.25)
  polygon(x=c(min.yr:max.yr,max.yr:min.yr),y=c(lower,upper[N.yr:1]),col="gray90",border="gray90")
  #points(years+(min.yr-1),weights,pch = ".")
  points(seq(min.yr,max.yr),fitted,type="l",col="red",lwd=2)
  
}
emf("SSB_StateSpace.emf",width = 14, height = 7,emfPlus = FALSE)
graph.ssm.SSB(exp(N.est),N.yr.ab,w,years.ab,max.yr.ab,min.yr.ab)
dev.off()

# lags in X (FCL) that predict Y (SSB):
ccf(fitted_FCL[1:46],fitted_SSB) # dang... that's rad.
ccfvalues<-ccf(fitted_FCL[1:46],fitted_SSB)
# Thus, FCL from 2 years ago best predicts SSB in current year (and it predicts it pretty darn well)
plot(fitted_FCL[1:46],exp(fitted_SSB))


###BRING IN RECRUITMENT DATA#######################################################################

REC<-read.csv("ancovy_recruitment.csv", header = T)
weights<-log(REC$Age.0.anchovy..10.3.mt.)

# Now I'm going to change the years so that the earliest year equals 1...
max.yr.r<-max(REC$Year); min.yr.r<-min(REC$Year)
years.r<-REC$Year-min(REC$Year) +1
Source<-as.numeric(REC$Source) #get sources (factors) as numric for jags (if needed)
# Get the length of the data set
N.ver.r <-dim(REC)[1]
#get the number of years the data set spans
N.yr.r<-max(years.r)
### LETS MODEL 


# Specify model in BUGS language ####################
sink("ssm.jags")
cat("
model { 
# Priors and constraints
N.est[1] ~ dunif(0, 4)            # Prior for initial yearly max size
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
jags.data <- list(weights=weights,years=years.r,N.ver=N.ver.r, N.yr=N.yr.r)

# Parameters monitored
parameters <- c("sigma2.obs", "sigma2.proc", "N.est")

# MCMC settings
ni <- 25000
nt <- 3
nb <- 10000
nc <- 3

# Call JAGS from R ###########
ssm.rec <- jags(jags.data, inits=NULL, parameters, "ssm.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
attach.jags(ssm.rec)
fitted_Rec<-apply(N.est,2,mean)
##

### Graph results ##############3
graph.ssm<-function(N.est,N.yr,weights,years,max.yr, min.yr){
  fitted<-lower<-upper<-numeric()
  
  for (i in 1:N.yr){
    fitted[i]<-mean(N.est[,i])
    lower[i]<-quantile(N.est[,i],0.27)
    upper[i]<-quantile(N.est[,i],0.75)
  }
  m1<-min(c(fitted,lower))
  m2<-max(c(fitted,upper)) 
  par(mar=c(4.5,4,1,1),cex=1.2)
  plot(0,0,ylim=c(m1,m2),xlim=c(1960,2005),ylab="Recruitment",xlab="year",las=1,col="black",type="l",lwd=2,frame=FALSE,axes=FALSE)
  axis(2,las=1)
  axis(1,at=seq(min.yr,max.yr,5), labels=seq(min.yr,max.yr,5))
  #axis(1,at=0:N.yr,labels=rep("",N.yr+1),tcl=-0.25)
  polygon(x=c(min.yr:max.yr,max.yr:min.yr),y=c(lower,upper[N.yr:1]),col="gray90",border="gray90")
  #points(years+(min.yr-1),weights,pch = ".")
  points(seq(min.yr,max.yr),fitted,type="l",col="red",lwd=2)
  
}
graph.ssm(exp(N.est),N.yr,weights,years.r,max.yr.r,min.yr.r)

ccf(fitted_FCL[5:46],fitted_Rec[1:42]) # specific cells to address year offsets of data (1964-2005)
##

###BRING IN LARVAE DATA#######################################################################

LAR<-read.csv("anchovy.larvae.by.year.csv", header = T)
weights<-log(LAR$ratioBig.Small)

# Now I'm going to change the years so that the earliest year equals 1...
max.yr.L<-max(LAR$Year); min.yr.L<-min(LAR$Year)
years.L<-LAR$Year-min(LAR$Year)+1

# Get the length of the data set
N.ver.L <-dim(LAR)[1]
#get the number of years the data set spans
N.yr.L<-max(years.L)
### LETS MODEL 


# Specify model in BUGS language ####################
sink("ssm.jags")
cat("
model { 
# Priors and constraints
N.est[1] ~ dunif(-10, 10)            # Prior for initial yearly max size
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
jags.data <- list(weights=weights,years=years.L,N.ver=N.ver.L, N.yr=N.yr.L)

# Parameters monitored
parameters <- c("sigma2.obs", "sigma2.proc", "N.est")

# MCMC settings
ni <- 25000
nt <- 3
nb <- 10000
nc <- 3

# Call JAGS from R ###########
ssm.lar <- jags(jags.data, inits=NULL, parameters, "ssm.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
attach.jags(ssm.lar)
fitted_Lar<-apply(N.est,2,median)
##

### Graph results ##############3
graph.ssm<-function(N.est,N.yr,weights,years,max.yr, min.yr){
  fitted<-lower<-upper<-numeric()
  
  for (i in 1:N.yr){
    fitted[i]<-median(N.est[,i])
    lower[i]<-quantile(N.est[,i],0.27)
    upper[i]<-quantile(N.est[,i],0.75)
  }
  m1<-min(c(fitted,lower))
  m2<-max(c(fitted,upper)) 
  par(mar=c(4.5,4,1,1),cex=1.2)
  plot(0,0,ylim=c(m1,m2),xlim=c(1960,2005),ylab="Ratio Big:Small",xlab="year",las=1,col="black",type="l",lwd=2,frame=FALSE,axes=FALSE)
  axis(2,las=1)
  axis(1,at=seq(min.yr,max.yr,5), labels=seq(min.yr,max.yr,5))
  #axis(1,at=0:N.yr,labels=rep("",N.yr+1),tcl=-0.25)
  polygon(x=c(min.yr:max.yr,max.yr:min.yr),y=c(lower,upper[N.yr:1]),col="gray90",border="gray90")
  #points(years+(min.yr-1),weights,pch = ".")
  points(seq(min.yr,max.yr),fitted,type="l",col="red",lwd=2)
  
}
emf("ratio_StateSpace.emf",width = 14, height = 7,emfPlus = FALSE)
graph.ssm(exp(N.est),N.yr.L,weights,years.L,max.yr.L,min.yr.L)
dev.off()

ccf(fitted_FCL,fitted_Lar[10:55]) # specific cells to address year offsets of data (1960-2009)

ccfvalues<-ccf(fitted_FCL,fitted_Lar[10:55])
##

