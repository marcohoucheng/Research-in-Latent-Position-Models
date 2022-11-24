rm(list=ls())
library(latentnet)
library(MASS)
library(STAT)
library(vegan)
library(randomcoloR)
data(sampson)
Names = samplike %v% "vertex.names"
samp.fit <- ergmm(samplike ~ euclidean(d=2))
#samp.fit <- ergmm(samplike ~ euclidean(d=2), control = control.ergmm(sample.size = 10000, burnin = 10000))

Monks = samp.fit$sample$Z
k = length(Monks[,1,1])
n = length(Monks[1,,1])
d = length(Monks[1,1,])

start_time <- Sys.time()
Monksrot = array(0, c(k,n,d))
for (i in 1:k){
  Monksrot[i,,] = procrustes(Monks[1,,], Monks[i,,], scale = F)$Yrot
}

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
AD = Y + t(Y)
CMD = cmdscale(AD, k = 2)
CMDpro = procrustes(OS, CMD, scale = F)
CMDstar = CMDpro$Yrot + cbind(rep(CMDpro$translation[,1],n),rep(CMDpro$translation[,2],n))
mkl = samp.fit$mkl$Z
pmode = samp.fit$mcmc.pmode$Z
pmean = array(0, c(n,d))
for (i in 1:n){
  for (t in 1:d){
    pmean[i,t] = mean(Monksrot[,i,t])
  }
}

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

OSdistances = distance(OS)
CMDdistances = distance(CMD)
mkldistances = distance(mkl)
pmodedistances = distance(pmode)
pmeandistances = distance(pmean)

OSscore = min(M)
CMDscore = scores(CMDdistances)
mklscore = scores(mkldistances)
pmodescore = scores(pmodedistances)
pmeanscore = scores(pmeandistances)

#Plot the scores
Scores = c(OSscore, CMDscore, mklscore, pmodescore, pmeanscore)
scorebar <- barplot(Scores,
        #main = "Scores",
        xlab = "Method",
        ylab = "Loss",
        ylim = c(0,max(Scores)+50),
        names.arg = c("OS", "CMD", "mkl", "pmode", "pmean"))
method <- c(OSscore, CMDscore, mklscore, pmodescore, pmeanscore)
text(scorebar,method+10,labels=as.character(round(method,2)))
grid()

#Simulation
#beta = array(0,c(n,n))
#betapmode = array(0,c(n,n))
#beta[lower.tri(beta)] <- samp.fit$mkl$beta
#betapmode[lower.tri(betapmode)] <- samp.fit$mcmc.pmode$beta
beta <- samp.fit$mkl$beta
betapmode <- samp.fit$mcmc.pmode$beta


probability <- function(b,distances){
  Exp = exp(b - distances)
  Exp[upper.tri(Exp, diag=T)] <- 0
  prob = Exp/(1+Exp)
  Prob = prob + t(prob)
  return(Prob)
}

simulation <- function(t,n,prob){
  totaledges = array(0,t)
  for (k in 1:t){
    unif = array(0,c(n,n))
    edges = array(0,c(n,n))
    for (i in 1:n){
      for (j in 1:n){
        unif[i,j]=runif(1)
      }
    }
    edges = ifelse((prob-unif) > 0, 1, 0)
    totaledges[k] = sum(edges)
  }
  return(mean(totaledges))
}

OSprob = probability(beta, OSdistances)
CMDprob = probability(beta, CMDdistances)
mklprob = probability(beta, mkldistances)
pmeanprob = probability(beta, pmeandistances)
pmodeprob = probability(betapmode, pmodedistances)

True_no_of_edges = sum(as.matrix(samplike))
OS_average_edge = simulation(10000,n,OSprob)
CMD_average_edge = simulation(10000,n,CMDprob)
mkl_average_edge = simulation(10000,n,mklprob)
pmode_average_edge = simulation(10000,n,pmodeprob)
pmean_average_edge = simulation(10000,n,pmeanprob)

True_no_of_edges
OS_average_edge
CMD_average_edge
mkl_average_edge
pmode_average_edge
pmean_average_edge

True_no_of_edges = sum(as.matrix(samplike))
OS_average_edge = sum(OSprob)
CMD_average_edge = sum(MDprob)
mkl_average_edge = sum(mklprob)
pmode_average_edge = sum(pmodeprob)
pmean_average_edge = sum(pmeanprob)

True_no_of_edges
OS_average_edge
CMD_average_edge
mkl_average_edge
pmode_average_edge
pmean_average_edge

end_time <- Sys.time()
end_time - start_time

#Plots
mklpro = procrustes(OS, mkl, scale = F)
mklstar = mklpro$Yrot + cbind(rep(mklpro$translation[,1],n),rep(mklpro$translation[,2],n))
pmodepro = procrustes(OS, pmode, scale = F)
pmodestar = pmodepro$Yrot + cbind(rep(pmodepro$translation[,1],n),rep(pmodepro$translation[,2],n))
pmeanpro = procrustes(OS, pmean, scale = F)
pmeanstar = pmeanpro$Yrot + cbind(rep(pmeanpro$translation[,1],n),rep(pmeanpro$translation[,2],n))
plot(OS, col = "white", xlab = "", ylab = "")
points(OS, pch = 20)
text(OS, labels = Names, pos = 4, cex = 0.5)
points(CMDstar, pch = 20, col = "red")
text(CMDstar, labels = Names, col = "red", pos = 4, cex = 0.5)
points(mklstar, pch = 20, col = "purple")
text(mklstar, labels = Names, col = "purple", pos = 4, cex = 0.5)
points(pmodestar, pch = 20, col = "blue")
text(pmodestar, labels = Names, col = "blue", pos = 4, cex = 0.5)
points(pmeanstar, pch = 20, col = "orange")
text(pmeanstar, labels = Names, col = "orange", pos = 4, cex = 0.5)

colour = distinctColorPalette(n)
plot(OS, col = colour, pch = 21, xlab = "", ylab = "")
text(OS, col = colour, labels = Names, pos = 4)
points(CMDstar, col = colour, pch = 22)
#text(CMDstar, col = colour, labels = Names, pos = 4, cex = 0.5)
points(mklstar, col = colour, pch = 23)
#text(mklstar, col = colour, labels = Names, pos = 4, cex = 0.5)
points(pmodestar, col = colour, pch = 24)
#text(pmodestar, col = colour, labels = Names, pos = 4, cex = 0.5)
points(pmeanstar, col = colour, pch = 25)
#text(pmeanstar, col = colour, labels = Names, pos = 4, cex = 0.5)

par(mfrow=c(3,2))
for (i in 1:6){
  plot(Monksrot[,i,], col = "grey", xlab="", ylab="", main = paste("Actor", i), pch = 20, cex = 0.1)
  points(OS[i,], col = "black")
  points(CMDstar[i,], col = "red")
  points(mklstar[i,], col = "purple")
  points(pmodestar[i,], col = "blue")
  points(pmeanstar[i,], col = "orange")
}

colour = distinctColorPalette(n)
par(mfrow=c(1,1))
plot(Monks[,1,], col=colour[i], xlab="", ylab="", main = "Actor 1", xlim = c(-5,5), ylim = c(-5,5), pch = 20, cex = 0.05)
for (i in 2:n){
  points(Monksrot[,i,], col = colour[i], xlab="", ylab="", main = paste("Actor", i), pch = 20, cex = 0.05)
}
points(OS,col="black",pch=4)
points(CMDstar, col = "red", pch = 9)
text(CMDstar, labels = 1:n, col="red")

start_time <- Sys.time()
end_time <- Sys.time()
end_time - start_time
