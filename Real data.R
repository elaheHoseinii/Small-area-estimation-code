###-------------------------------------------------------------------

# R code-Intrinsic CAR-semiparametric POISS version NONE Spline by Real Data
rm(list=ls())

path<-getwd()

#install.packages("INLA", repos="https://inla.r-inla-download.org/R/stable")
library (INLA)
library(lme4)
library(igraph)
library(psych)
library(rgdal)
library(raster)
library(spdep)
library(igraph)
library(mclcar)
library(mvtnorm)
library(coda)
library(BRugs)
library(rjags)
library(R2OpenBUGS)
library(dclone)
library(stats)
library(nlme)
library(MASS)
library(splines)
library(lattice)
#library(ggplot2) # plots results
#library(RColorBrewer)  # for plot colour palette
#library(tmap)
#library(spatialreg)

#Data (insurance premium in IRAN):

#observed:
y1 = read.csv("C:/Users/Alireza/Desktop/Miss Hosseini/data2/y.csv", header = F)
X = read.csv("C:/Users/Alireza/Desktop/Miss Hosseini/data2/x.csv", header = F)
index = read.csv("C:/Users/Alireza/Desktop/Miss Hosseini/data2/index.csv", header = F)
index=index[,1]
y1=y1[,9]
X=X[,9]

plot(index, y1)
length(y1)
plot(X, y1)
###########################################################
#x <- bs(X, knots=quantile(X, c(0.1,.8,.9)), degree=3) # By Spline                                                 # Simple Linear Regression
#mm<-lm(y1~x)
#pred<-predict(mm)
#plot(sort(y1), type = "l", lwd=2, ylab = "Fit", xlab = expression(x))
#lines(sort(pred), col = 2, lwd=2, lty=2)
############################################################
#Population of Provinces in Iran:
E = read.csv("C:/Users/Alireza/Desktop/Miss Hosseini/data2/E.csv", header = F)
E= E[,1]
n=31

####### a spatial generalized semiparametric model with a ICAR spatial effect

###------------------------- Essential function -----------------------###
### define space domain based on Minnesota shape file
G.shape = readOGR(dsn=path.expand("C:\\Users\\Alireza\\Desktop\\Miss Hosseini\\iran-administrative\\IRAN PROVINC\\IRAN PROVINC\\ostan_1392.shp"),layer="ostan_1392")
names(G.shape)

plot(G.shape)

nb_q <- poly2nb(G.shape)

ccMat=nb2mat(nb_q,style='B') ### weight matrix
#is.symmetric.matrix(ccMat)
num = as.vector(rowSums(ccMat))


##---- define adjacency matrix to use in OpenBUGS
n.site <- dim(ccMat)[1]          # n.site: number of areas
n.edge <- sum(ccMat)/2		   # n.edge: number of unique pairs

# Find neighbouring areas

SEind1 <- SEind2 <- 0
matmy <- ccMat
for(i in 1:(n.site-1)){
  for(j in (i+1):n.site){
    if (ccMat[i,j]>0) {SEind1<-c(SEind1,i)
    SEind2<-c(SEind2,j)
    matmy[i,j]<-matmy[j,i]<-length(SEind1)-1
    }
  }
}

SEind1 <- SEind1[-1] # edges sorted by row of upper triangle of the adj. matrix
SEind2 <- SEind2[-1] # SEind1[k]=i and SEind2[k]=j => kth edge is edge ij

#dput(SEind1,"SEind1.txt")
#dput(SEind2,"SEind2.txt")
#dput(ccMat, "W.txt")

# create adjacency information needed for WinBUGS
mkAdj <- function(W){
  n <- nrow(W)
  adj <- 0
  for(i in 1:n){
    for(j in 1:n){
      if(W[i,j]==1){adj<-append(adj,j)
      }
    }
  }
  adj <- adj[-1]
  return(adj)
}

adj = mkAdj(ccMat)

###--------------------------- Simulation
set.seed(123456)
M = diag(num)
C = ccMat/num
n = nrow(C)

rho = 1
x = X
m = log(E) # log E_i
data.sim <- function(n=31, tau.phi){
  Phi <- rnorm(n)
  Phi.Mean <- c()
  for (i in 1:n){
    Phi.Mean[i] <- rho*(sum(C[i,-i]*Phi[-i])/num[i])
  }
  Phi.Sig <- 1/(tau.phi*num)
  Phi <- rnorm(n, Phi.Mean, Phi.Sig)
 
  eta = m + Phi
  lambda = exp(eta)
  y = sample(y1, length(y1), replace = T) # y1 is real data
return(list(O=y, eta = eta ,E=E, lambda = lambda, u = Phi))
}



rep=1
y.rep = eta.rep = lambda.rep = u = l.rep = eta.hat.ebrep = matrix (NA, ncol = rep, nrow = n)
############################## REP for all####################################################
for (q in 1:rep) {
  dt <- data.sim(n=31, tau.phi = 0.5)
  y.rep[,q] = dt$O
  eta.rep[,q] = dt$eta
  lambda.rep[,q] = dt$lambda
  u[,q] = dt$u

################ Fitting SGLMM with proper CAR using HB method
model <- function(){
# =====
# Likelihood
for(k in 1:K){
for (i in 1 : N) {
  O[i,k]  ~ dpois(mu[i,k])
  log(mu[i,k]) <- log(E[i]) + inprod(X[i,],beta[1,]) + S[k,i]
# Area-specific relative risk 
# RR[k,i] <- exp(S[k,i])
  }
}

# Proper CAR prior distribution for spatial random effects: 
for(k in 1:K){
 S[k,1:N] ~ car.normal(adj[], weights[], num[], prec)
}
  for(l in 1:sumNumNeigh) {
	weights[l] <- 1
}			

# Other priors:
# alpha  ~ dnorm(0, 0.0001)  
for (j in 1:np) {
 beta[1,j] ~ dnorm(0, 0.0001)
}
# prior on precision
prec  ~ dgamma(0.05, 0.05)
sigma <- sqrt(1/prec)		# standard deviation
#S.mean <- sum(S[])
}

#=============================================
#Data

N = length(y.rep[,q])
O = y.rep[,q]
E = dt$E

X <- bs(m, knots=quantile(m, c(0.25,.5,.75)), degree=3) # By Spline
#X <- m                                                 # Simple Linear Regression
X <- model.matrix(~X)
np = ncol(X)
num = num 
sumNumNeigh = sum(num)
######---------------

#adjacency matrix
adj = adj


K <- 1
Y <- matrix(c(O), N, K, byrow=F)
#parameters.hb <- c("beta", "prec", "sigma")
#data.hb <- list(O=Y, X=as.matrix(X), E=E, N=N, adj=adj, num=num, 
  #              sumNumNeigh=sumNumNeigh, K=1, np=np)

#initial values:
#inits.hb <- list(list(beta=matrix(rep(0, np),1,np), prec=1),
 #                list(beta=matrix(rep(0, np),1,np), prec=1))

#modelout.hb <- bugs.fit(data=data.hb, parameters.hb, model,
   #                     inits=inits.hb, n.chains=2, n.iter=3000, n.thin=5, n.burnin=1000, 
   #                     program="openbugs", debug=FALSE)

#print(paste("Gelman"))
#print(max((gelman.diag(modelout.hb))$psrf[,1]))
#stats.hb <- summary(modelout.hb,quantiles=c(0.005,0.01,0.025,0.05,0.5,0.95,0.975,0.99,0.995))

#lambda.hb <- as.numeric(lambdamax.diag(modelout.hb))
#prec.hb = stats.hb$statistics[np+1,1]
#sigma.hb = stats.hb$statistics[np+2,1]
#Beta.hb = stats.hb$statistics[1:np,1]
#==================================================

#DC approach:
dat2 <- list(O=dcdim(O), X=X, E=E, N=N, adj=adj, num=num, 
             sumNumNeigh=sumNumNeigh, np=np, K=1)
parameters.dc <- c("beta", "prec", "sigma")
inits.dc <- list(list(beta=matrix(rep(0, np),1,np), prec=1),
                 list(beta=matrix(rep(0, np),1,np), prec=1))

K1 <- c(0) #number of clones
dcmod <-dc.fit(dat2, parameters.dc, model, inits.dc, n.iter=5000, 
               n.burnin=2000, n.thin=10, n.chains=2, n.clones=c(3), multiply=c("K"),	
               unchanged=c("N", "X", "np", "sumNumNeigh", "E", "num", "adj"), 
               flavour="bugs", program="openbugs", debug=TRUE)

dcd <- dcdiag(dcmod)
r2 <- dcd[1,4]
mse <- dcd[1,3]
r.hat <- dcd[1,5]
#lambda.ratio <- dcd[1,2]/lambda.hb
stats.dc = summary(dcmod)
sigma.dc = stats.dc$statistics[np+2,1]
prec.dc = stats.dc$statistics[np+1,1]
Beta.dc = stats.dc$statistics[1:np,1]

#=====================prediction part is deleted==================================================
X.Beta = X%*%Beta.dc
##########################################################
eta0=c()
for(i in 1: nrow(l.rep)){
  eta0[i]<-eta.rep[which.max(y.rep[i,])]
}

sigma = (sigma.dc * num)*(diag(rep(1,31))) # Sigma Matrix
P<- matrix(0, ncol=nrow(eta.rep), nrow= nrow(eta.rep))
dim(P)
for(j in 1:ncol(P)){
  for(i in 1:nrow(P)){
    if(i==j){ P[i,j]<-(1/exp(eta0[i]))}
  }
  P
}

Z<- diag(ncol(sigma))

R<- Z%*%sigma%*%t(Z)+P
#..........................................................

######### tolid eta hat EB marhale1
#for(q in 1: ncol(l.rep)){
for(i in 1: nrow(l.rep)){
  l.rep[i,q]<- (y.rep[i,q] - exp(eta0[i])+eta0[i]%*%exp(eta0[i]))/exp(eta0[i])
  #  }
}
l.rep
#==
dim(y.rep)
#for(q in 1: ncol(eta.hat.ebrep)){
for(i in 1:nrow(eta.hat.ebrep)){
  eta.hat.ebrep[i,q]<- X.Beta[i]+ t(Z[i,])%*%sigma%*%t(Z)%*%solve(R)%*%(l.rep[,q]-X.Beta)
  #}
}
eta.hat.ebrep
}
########################################################
#=======================================================
########################################################

mspe<- function(etahat, etahat.eb){
  mspe<- c(rep(0, nrow(u)))
  for(i in 1:nrow(etahat)){
    for(j in 1: ncol(etahat))
      mspe[i]<-sum(etahat[i,j]-etahat.eb[i,j])^2/ncol(etahat)
  }
  return(mspe)
}

MSE.Spline<- mspe(eta.rep,eta.hat.ebrep)
boxplot(MSE.Spline)

save(MSE.Spline, file = "MSE-Spline.RData")
save(eta.rep, file = "eta-Spline.RData")
save(eta.hat.ebrep, file = "eta-Spline.RData")



