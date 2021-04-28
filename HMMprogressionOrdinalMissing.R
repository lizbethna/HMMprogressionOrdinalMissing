##################################################
### Paper:
### Monitoring Parkinson’s disease progression based on recorded speech 
### with missing ordinal responses and replicated covariates
###
### Authors:
### Lizbeth Naranjo (1), Carlos J. Perez (2), Yolanda Campos-Roca (3).
###
### (1) Departamento de Matemáticas, Facultad de Ciencias, 
### Universidad Nacional Autonoma de Mexico (UNAM), Mexico
### (2) Departamento de Matematicas, Facultad de Veterinaria, 
### Universidad de Extremadura, Spain
### (3) Departamento de Tecnologias de los Computadores y de las Comunicaciones, 
### Escuela Politecnica, Universidad de Extremadura, Spain
### 
### Journal: 
### Submitted. Under Revision. 
##################################################

##################################################
### R packages
###
### Load the R libraries
##################################################
library(rjags)
#library(mnormt)
library(MCMCpack) ### MCMC
#library(LearnBayes) ### Bayes
#library(gamlss.dist) ### Inverse Gaussian distribution
#library(mnormt)
library(coda)
#library(dclone) # To run MCMC in
#library(snow)   # several cores
#library(corrplot)
#library(pROC)
#library(xlsx)
#library(HDInterval)

##################################################

##################################################
### ADDRESS
###
### Instructions: 
### Change the address where the data and codes are located. 
### setwd("HERE")
##################################################
setwd("~/Documents/Articulos/ReplicacionesORL/ArticuloProgresionHY/HMMordinalMissingCodes/") 
getwd()
##################################################

##################################################
### SIMULATE DATA
###
### Run the followin lines of code.
##################################################

##################################################
### Truncated Normal I[tra < Z]
rnormleft <- function(tra,mu,sig){
  rp <- pnorm(tra, mean=mu, sd=sig)
  u <- rp+(1-rp)*runif(1)
  q <- qnorm(u, mean=mu, sd=sig)
  if(!is.finite(q)){ q = tra }
  return(q)
}
##################################################

##################################################
### Parameters 

N = 100   # subjects
TT = 4   # times
K = 4   # categories
J =  3   # replicaciones
L = 4   # covariates with replications
M = 2   # covariates without replications

beta = c(1.0,0.5,-1.5,-2)
gamma = c(0.4,0.2)
kappa = c(-Inf,0,1.5,3.0,Inf)
sigma2 = c(0.2,0.25,0.3,0.35)^2

##################################################

##################################################
### Covariates and Replications 

Z = array( NA ,dim=c(N,M)) 
Z[,1] = rbinom(N,prob=0.5,size=1)
Z[,2] = rpois(N,lambda=5)

Xmean = array( runif(N*TT*L) ,dim=c(N,TT,L)) 
X = array(NA,dim=c(N,TT,J,L))
for(i in 1:N){
for(t in 1:TT){	
for(l in 1:L){
	X[i,t,,l] = rnorm(J,Xmean[i,t,l],sqrt(sigma2[l]))
}	}	}

##################################################

##################################################
### Simulate Ordinal Responses

set.seed(12345)

eta = array(NA,dim=c(N,TT))   ### linear predictor
eta[,1] = Xmean[,1,]%*%beta + Z[,]%*%gamma 
for(t in 2:TT){
  eta[,t] = Xmean[,t,]%*%beta + Z[,]%*%gamma 
}
Y = W = array(NA,dim=c(N,TT))
for(i in 1:N){
  W[i,1] = rnorm(1,eta[i,1],1)
  for(t in 2:TT){
    W[i,t] = rnormleft(W[i,t-1],eta[i,t],1)
}  }
for(i in 1:N){
  for(k in 1:K){
    if(kappa[k]< W[i,1] & W[i,1]<=kappa[k+1]){
      Y[i,1] = k
  }  }
  for(t in 2:TT){
    for(k in 1:K){
      if(kappa[k]< W[i,t] & W[i,t]<=kappa[k+1]){
        Y[i,t] = k
}  }  }  }	
  
Ytrue = Y
Wtrue = W
Xtrue = X
  
plot(c(1:K),W[1,], main="W", type="n",col=1, ylim=c(min(W),max(W)))
for(i in 1:N){	lines(c(1:K),W[i,],col=i)		}
for(t in 2:TT){ print(table(Y[,t-1],Y[,t])) }
tabla = matrix(0,K,K)
for(i in 1:N){
  for(t in 2:TT){ 
    tabla[Y[i,t-1],Y[i,t]] = tabla[Y[i,t-1],Y[i,t]] +1
}  }
tabla

##################################################

##################################################
### Simulate the Missing Data

NTmiss = N*TT*0.10   ### number of missing data
missN = sample(c(1:N),size=NTmiss,replace=TRUE)
missT = sample(c(2:TT),size=NTmiss,replace=TRUE)
for(i in 1:NTmiss){
  Y[missN[i],missT[i]] = NA
  X[missN[i],missT[i],,] = NA
}

##################################################

##################################################
### Identify the missing NA's 

Ni = rep(NA,N)   ### Indica cuantos valores son distintos de NA, para cada sujeto i=1:N
Nx = Nx.noNA = Nx.NA = array(NA,dim=c(N,TT))   

for(i in 1:N){
  Ni[i] = sum(!is.na(Y[i,1:TT]))  ### denotes how many values are not NA's
  for(t in 1:TT){
    Nx[i,t] = ifelse(!is.na(Y[i,t]),t,NA)   ### observed years
    Nx.NA[i,t] = ifelse(is.na(Y[i,t]),t,NA)   ### denotes which years for subject are NA's, t=1:Ni, i=1:N.
  }
  Nx.noNA[i,1:Ni[i]] = Nx[i,!is.na(Nx[i,])]  ### denotes which years ser not NA's, t=1:Ni, i=1:N.
}  


N.mis = sum(Ni<TT)
id.mis = Ni.mis = rep(NA,N.mis)
TTi.mis = matrix(NA,N.mis,TT)
id = 0
for(i in 1:N){
  if(Ni[i]<TT){
    id = id+1
    id.mis[id] = i
    Ni.mis[id] = TT-Ni[i]
    TTi.mis[id,1:Ni.mis[id]] = Nx.NA[i,!is.na(Nx.NA[i,])]
  }
}  

##################################################

##################################################
### DATA, PARAMETERS, INITIAL VALUES

### Data     
data <- list(
  Y = Y-1 , ### caterogies are 0,1,2,...,K-1
  X = X , 
  Z = Z ,
  N = N , 
  Ni = Ni ,
  J = J ,
  K = K ,
  L = L ,
  M = M ,
  TT = TT ,
  Nx.noNA = Nx.noNA ,
  N.mis = N.mis ,
  id.mis = id.mis , 
  Ni.mis = Ni.mis ,
  TTi.mis = TTi.mis  
)

### Initial Values
Wobs =  cbind(rep(-1.5,N),Y-1.5,rep(6.5,N))
for(i in 1:N){
  for(t in 1:TT){
    if(is.na(Wobs[i,t+1]))
    Wobs[i,t+1] = mean(c(Wobs[i,t],Wobs[i,t+2]),na.rm=TRUE)
  }
}
Wobs[,1] = NA
Wobs[,TT+2] = NA
inits <- function(){	list(
  "beta" = rnorm(L,0,0.001) , 
  "gamma" = rnorm(M,0,0.001) ,
  "invsigma2" = rep(1,L) ,
  "kappaini" = c(0:(K-2)) , 
  "Wobs" = Wobs 
)	}

### Parameters
param <- c(
  "beta","gamma", "sigma2" , "kappa"
)

##################################################

##################################################
### FIT THE MODEL

### Fit the Model
fit.model <- jags.model("HMMprogressionOrdinalMissing.bug", data, inits, n.chains=3)

update(fit.model,10000)

sample.ord <- coda.samples(fit.model, param,  n.iter=10000, thin=5)

### Summary & Graphics
plot(sample.ord)
(post.ord <- summary(sample.ord))

### Checking convergence quickly
#print(gelman.diag(sample.ord[,param[]]))
print(100*batchSE(sample.ord, batchSize=50)/post.ord$statistics[,2])
print(effectiveSize(sample.ord))


##################################################
##################################################
### Graphics


kappaF = c(post.ord$statistics["kappaF[1]",1], post.ord$statistics["kappaF[2]",1], post.ord$statistics["kappaF[3]",1])
kappaP = c(post.ord[1]$statistics["kappaP[1]",1], post.ord$statistics["kappaP[2]",1], post.ord$statistics["kappaP[3]",1])
YfitF = c(kappaF[1]-(kappaF[2]-kappaF[1])/2,(kappaF[1]+kappaF[2])/2,(kappaF[2]+kappaF[3])/2,kappaF[3]+(kappaF[3]-kappaF[3])/2)
YfitP = c(kappaP[1]-(kappaP[2]-kappaP[1])/2,(kappaP[1]+kappaP[2])/2,(kappaP[2]+kappaP[3])/2,kappaP[3]+(kappaP[3]-kappaP[2])/2)


par(mfrow=c(3,3))
for(i in 1:18){
  plot(c(1,2,3,4),
       c(post.ord$statistics[paste0("Wobs[",i,",2]"),1], post.ord$statistics[paste0("Wobs[",i,",3]"),1], post.ord$statistics[paste0("Wobs[",i,",4]"),1], post.ord$statistics[paste0("Wobs[",i,",5]"),1]),
       main=paste("Subject ",i), type="l",col="black",lwd=2,ylim=c(0-max(kappaF[2],kappaP[2]),max(kappaF[3],kappaP[3])+max(kappaF[3]-kappaF[2],kappaP[3]-kappaP[2])),
       xlim=c(1,4),
       xlab="Time",ylab="Stages")
  lines(c(1,2,3,4),
        c(post.ord[2]$quantiles[paste0("Wobs[",i,",2]"),1], post.ord$quantiles[paste0("Wobs[",i,",3]"),1], post.ord$quantiles[paste0("Wobs[",i,",4]"),1], post.ord$quantiles[paste0("Wobs[",i,",5]"),1]),
        col="black",lty=2,lwd=2)
  lines(c(1,2,3,4),
        c(post.ord$quantiles[paste0("Wobs[",i,",2]"),5], post.ord$quantiles[paste0("Wobs[",i,",3]"),5], post.ord$quantiles[paste0("Wobs[",i,",4]"),5], post.ord$quantiles[paste0("Wobs[",i,",5]"),5]),
        col="black",lty=2,lwd=2)
  points(c(1,2,3,4),
         c(YfitF[Y[i,1]],YfitP[Y[i,2]],YfitP[Y[i,3]],YfitP[Y[i,4]]),
         cex=1.5,col="blue",pch=19)
  abline(h=kappaP,lty=3,col="gray" )
}


##################################################
##################################################

##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
