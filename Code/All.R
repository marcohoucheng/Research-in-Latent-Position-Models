rm(list=ls())
library(randomcoloR)
library(STAT) #cmdscale
library(vegan)
library(latentnet)

n = 10 #no. of actors
d = 2 #latent spaces' dimensions
mu = 0 #mean
var = 1 #small variance
k = 1000 #number of iterations
Var = 10 #big variance

GaussianLatentSpaceGenerator = function(n,d,mu,var,k,Var){ 
  True = replicate(d, rnorm(n, mu, Var))
  S = replicate(k, True)
  for (i in 1:k){
    e = replicate(d, rnorm(n, mu, var))
    S[,,i] = True + e
  }
  colour = distinctColorPalette(n)
  plot(S[,,1], col=colour, xlab="", ylab="", xlim = c(min(S[,1,])-0.5*Var,max(S[,1,])+0.5*Var), ylim = c(min(S[,2,])-0.5*Var,max(S[,2,])+0.5*Var), pch = 20)
  for (i in 2:k){
    points(S[,,i],col=colour, pch = 20)
  }
  #text(S[,,1], labels = 1:n)
  assign("S", S, envir=globalenv())
  assign("True", True, envir=globalenv())
}

GaussianLatentSpaceGenerator(n,d,mu,var,k,Var) #Generate Positions

#Here I'm recording the time and the accuracy for the 3 methods, one by one
start_time <- Sys.time()
M = array(0, k)
X = array(0, c(n,n))
Y = array(0, c(n,n))
Z = array(0, c(n,n,k))

for (i in 1:n){
  for (j in 1:i){
    for (t in 1:k){
      Z[i,j,t] = sqrt(sum((S[i,,t]-S[j,,t])^2))
    }
    Y[i,j] = mean(Z[i,j,])
  }
}

for (m in 1:k){
  for (i in 1:n){
    for (j in 1:i){
      X[i,j] = (Z[i,j,m]-Y[i,j])^2
    }
  }
  M[m] = sum(X)
}
score1 = min(M) #First accuracy

end_time <- Sys.time()
time1 = end_time - start_time #First method time taken

start_time <- Sys.time()

Y = array(0, c(n,n))
Z = array(0, c(n,n,k))

for (i in 1:n){
  for (j in 1:i){
    for (t in 1:k){
      Z[i,j,t] = sqrt(sum((S[i,,t]-S[j,,t])^2))
    }
    Y[i,j] = mean(Z[i,j,])
  }
}

AD = array(0, c(n,n))
AD = Y + t(Y)
CMD = cmdscale(AD, k = 2)

C = array(0, c(n,n))
D = array(0, c(n,n))
for (i in 1:n){
  for (j in 1:i){
    C[i,j] = sqrt(sum((CMD[i,]-CMD[j,])^2))
    D[i,j] = (C[i,j]-Y[i,j])^2
  }
}
score2 = sum(D) #Second accuracy

end_time <- Sys.time()
time2 = end_time - start_time #Second method time taken

start_time <- Sys.time()

Y = array(0, c(n,n))
Z = array(0, c(n,n,k))

for (i in 1:n){
  for (j in 1:i){
    for (t in 1:k){
      Z[i,j,t] = sqrt(sum((S[i,,t]-S[j,,t])^2))
    }
    Y[i,j] = mean(Z[i,j,])
  }
}

Srot = array(0, c(n,d,k))
for (i in 1:k){
  rot = procrustes(S[,,1], S[,,i], scale = F)
  Srot[,,i] = rot$Yrot + cbind(rep(rot$translation[,1],n),rep(rot$translation[,2],n))
}

PosMean = array(0, c(n,d))
for (i in 1:n){
  for (t in 1:d){
    PosMean[i,t] = mean(Srot[i,t,])
  }
}
E = array(0, c(n,n))
G = array(0, c(n,n))
for (i in 1:n){
  for (j in 1:i){
    E[i,j] = sqrt(sum((PosMean[i,]-PosMean[j,])^2))
    G[i,j] = (E[i,j]-Y[i,j])^2
  }
}
score3 = sum(G) #Third accuracy

end_time <- Sys.time()
time3 = end_time - start_time #Third method time taken

#Plots
colour = distinctColorPalette(n)
plot(S[,,1], col=colour, xlab="", ylab="", xlim = c(min(S[,1,])-0.5*Var,max(S[,1,])+0.5*Var), ylim = c(min(S[,2,])-0.5*Var,max(S[,2,])+0.5*Var), pch = 20, xaxt='n', yaxt='n',ann=FALSE)
for (i in 2:k){
  points(S[,,i],col=colour, pch = 20)
}
points(S[,,which.min(M)],col="black",pch = 20, cex = 1.5)
legend("topright", inset = 0.02, legend = "Bayes Estimator", pch = 20, cex = 2)
CMDpro = procrustes(True, CMD, scale = F)
CMDstar = CMDpro$Yrot + cbind(rep(CMDpro$translation[,1],n),rep(CMDpro$translation[,2],n))
points(CMDstar, pch = 17, cex = 1.5)
legend("topright", inset = 0.02, legend = c("Bayes Estimator", "CMD"), pch = c(20, 17), cex = 2)
PosMeanpro = procrustes(True, PosMean, scale = F)
PosMeanstar = PosMeanpro$Yrot + cbind(rep(PosMeanpro$translation[,1],n),rep(PosMeanpro$translation[,2],n))
points(PosMeanstar, pch = 15, cex = 1.5)
legend("topright", inset = 0.02, legend = c("Bayes Estimator", "CMD", "Average Position"), pch = c(20, 17, 15), cex = 2)

scores = c(score1, score2, score3)
times = c(0.1928961, 0.1603458, 0.650269)

par(mfrow=c(1,2))
scorebar <- barplot(times, names.arg = c("Bayes Estimator", "CMD", "Average Position"), main = "Loading time (seconds)", ylim = c(0,0.7), cex.names = 1.5, cex.axis = 1.5, cex.main = 1.5)
text(scorebar,times+0.02,labels=as.character(round(times,2)), cex = 1.5)
scorebar <- barplot(scores, names.arg = c("Bayes Estimator", "CMD", "Average Position"), main = "Loss", ylim = c(0,20), cex.names = 1.5, cex.axis = 1.5, cex.main = 1.5)
text(scorebar,scores+0.4,labels=as.character(round(scores,2)), cex = 1.5)

#Monks ##################################################################################################

rm(list=ls())
data(sampson)
Names = samplike %v% "vertex.names"
samp.fit <- ergmm(samplike ~ euclidean(d=2))
#samp.fit <- ergmm(samplike ~ euclidean(d=2), control = control.ergmm(sample.size = 10000, burnin = 10000))

Monks = samp.fit$sample$Z
k = length(Monks[,1,1])
n = length(Monks[1,,1])
d = length(Monks[1,1,])

#Score function
#Scores
distance <- function(Positions){
  dist = array(0, c(n,n))
  for (i in 1:n){
    for (j in 1:i){
      dist[i,j] = sqrt(sum((Positions[i,]-Positions[j,])^2))
    }
  }
  return(dist)
}
scores <- function(distances){
  return(sum((distances-Y)^2))
}

#Bayes Estimator
start_time <- Sys.time()
M = array(0, k)
X = array(0, c(n,n))
Y = array(0, c(n,n))
Z = array(0, c(n,n,k))

for (i in 1:n){
  for (j in 1:i){
    for (t in 1:k){
      Z[i,j,t] = sqrt(sum((Monks[t,i,]-Monks[t,j,])^2))
    }
    Y[i,j] = mean(Z[i,j,])
  }
}
for (m in 1:k){
  for (i in 1:n){
    for (j in 1:i){
      X[i,j] = (Z[i,j,m]-Y[i,j])^2
    }
  }
  M[m] = sum(X)
}

OS = Monks[which.min(M),,]
OSdistances = distance(OS)
OSscore = min(M)

end_time <- Sys.time()
time1 = end_time - start_time

#CMD
start_time <- Sys.time()
Y = array(0, c(n,n))
Z = array(0, c(n,n,k))

for (i in 1:n){
  for (j in 1:i){
    for (t in 1:k){
      Z[i,j,t] = sqrt(sum((Monks[t,i,]-Monks[t,j,])^2))
    }
    Y[i,j] = mean(Z[i,j,])
  }
}
AD = Y + t(Y)
CMD = cmdscale(AD, k = 2)
CMDpro = procrustes(OS, CMD, scale = F)
CMDstar = CMDpro$Yrot + cbind(rep(CMDpro$translation[,1],n),rep(CMDpro$translation[,2],n))
CMDdistances = distance(CMD)
CMDscore = scores(CMDdistances)

end_time <- Sys.time()
time2 = end_time - start_time

#Position Mean
start_time <- Sys.time()
Y = array(0, c(n,n))
Z = array(0, c(n,n,k))

for (i in 1:n){
  for (j in 1:i){
    for (t in 1:k){
      Z[i,j,t] = sqrt(sum((Monks[t,i,]-Monks[t,j,])^2))
    }
    Y[i,j] = mean(Z[i,j,])
  }
}
Monksrot = array(0, c(k,n,d))
for (i in 1:k){
  Monksrot[i,,] = procrustes(Monks[1,,], Monks[i,,], scale = F)$Yrot
}
pmean = array(0, c(n,d))
for (i in 1:n){
  for (t in 1:d){
    pmean[i,t] = mean(Monksrot[,i,t])
  }
}

pmeandistances = distance(pmean)
pmeanscore = scores(pmeandistances)
end_time <- Sys.time()
time3 = end_time - start_time

#The other 2
mkl = samp.fit$mkl$Z
mkldistances = distance(mkl)
mklscore = scores(mkldistances)
pmode = samp.fit$mcmc.pmode$Z
pmodedistances = distance(pmode)
pmodescore = scores(pmodedistances)

#Plot the scores
par(mfrow=c(1,2))
times = c(1.466274,1.297034,2.709907)
scorebar <- barplot(times, names.arg = c("Bayes Estimator", "CMD", "Average Position"), main = "Loading time (seconds)", ylim = c(0,3), cex.names = 1.5, cex.axis = 1.5, cex.main = 1.5)
text(scorebar,times+0.075,labels=as.character(round(times,2)), cex = 1.5)
score = c(OSscore, CMDscore, pmeanscore, mklscore, pmodescore)
scorebar <- barplot(score, names.arg = c("Bayes Estimator", "CMD", "Average Position", "MKL", "Posterior Mode"), main = "Loss", ylim = c(0,450), cex.names = 1.2, cex.axis = 1.5, cex.main = 1.5)
text(scorebar,score+10,labels=as.character(round(score,2)), cex = 1.5)

#Simulation
beta <- samp.fit$mkl$beta
betapmode <- samp.fit$mcmc.pmode$beta
probability <- function(b,distances){
  Exp = exp(b - distances)
  Exp[upper.tri(Exp, diag=T)] <- 0
  prob = Exp/(1+Exp)
  Prob = prob + t(prob)
  return(Prob)
}

OSprob = probability(beta, OSdistances)
CMDprob = probability(beta, CMDdistances)
mklprob = probability(beta, mkldistances)
pmeanprob = probability(beta, pmeandistances)
pmodeprob = probability(betapmode, pmodedistances)
True_no_of_edges = sum(as.matrix(samplike))
OS_average_edge = sum(OSprob)
CMD_average_edge = sum(CMDprob)
mkl_average_edge = sum(mklprob)
pmode_average_edge = sum(pmodeprob)
pmean_average_edge = sum(pmeanprob)
average_edges = c(OS_average_edge, CMD_average_edge, pmean_average_edge, mkl_average_edge, pmode_average_edge)

par(mfrow=c(1,1))
scorebar <- barplot(average_edges, names.arg = c("Bayes Estimator", "CMD", "Average Position", "MKL", "Posterior Mode"), main = "Expected Number of Total Edges", ylim = c(0,100))
text(scorebar,average_edges+2,labels=as.character(round(average_edges)))
abline(h=True_no_of_edges, col = "Red", lty = 5, lwd = 2.5)
