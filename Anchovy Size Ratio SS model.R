 rm(list=ls()) #do you want to nuke the environ?
library(R2jags)
library(coda)
 library(devEMF)


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
graph.ssm(exp(N.est),N.yr.L,weights,years.L,max.yr.L,min.yr.L)

emf("ratio_StateSpace.emf",width = 14, height = 7,emfPlus = FALSE)
graph.ssm(exp(N.est),N.yr.L,weights,years.L,max.yr.L,min.yr.L)
dev.off()

# Now write results to file
fitted<-lower95<-upper95<-numeric()
for (i in 1:N.yr.L){
  fitted[i]<-mean(N.est[,i])
  lower95[i]<-quantile(N.est[,i],0.025)
  upper95[i]<-quantile(N.est[,i],0.975)
}
Year<-seq(from = min.yr.L, to = max.yr.L, by = 1) 
ratio_SS_MAP<- fitted
ratio_SS_lower95<- lower95
ratio_SS_upper95<- upper95

ratio.SS.estimates<-cbind(Year,ratio_SS_MAP,ratio_SS_lower95,ratio_SS_upper95)
write.csv(ratio.SS.estimates, file = "ratio.SS.estimates.csv")

ccf(fitted_FCL,fitted_Lar[10:55]) # specific cells to address year offsets of data (1960-2009)

ccfvalues<-ccf(fitted_FCL,fitted_Lar[10:55])
##

