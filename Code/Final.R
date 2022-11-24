rm(list=ls())
load("~/Desktop/GaussianLatentSpace.RData")
library(randomcoloR)
library(expm)
library(MASS) #sammon
library(STAT) #cmdscale
library(vegan)

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

GaussianLatentSpaceGenerator(n,d,mu,var,k,Var)

#Loss function 1 compare actor distances between samples

ML = array(0, k)
XL = array(0, c(n,n))
YL = array(0, c(n,n))
ZL = array(0, c(n,n,k))

for (i in 1:n){
  for (j in 1:i){
    for (t in 1:k){
      ZL[i,j,t] = sum((S[i,,t]-S[j,,t])^2)
    }
    YL[i,j] = sum(sqrt(ZL[i,j,]))
  }
}

for (m in 1:k){
  for (i in 1:n){
    for (j in 1:i){
      XL[i,j] = k*ZL[i,j,m]-2*sqrt(ZL[i,j,m])*YL[i,j]
    }
  }
  ML[m] = sum(XL)
}

#Loss function 2 compare actor distances to the mean distances

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

#Construct the "optimal" latent positions from the average distances

AD = array(0, c(n,n))
AD = Y + t(Y)
CMD = cmdscale(AD, k = 2)


A = array(0, c(n,n))
B = array(0, c(n,n))

for (i in 1:n){
  for (j in 1:i){
    A[i,j] = sum((CMD[i,]-CMD[j,])^2)
    B[i,j] = k*A[i,j]-2*sqrt(A[i,j])*YL[i,j]
  }
}

C = array(0, c(n,n))
D = array(0, c(n,n))
for (i in 1:n){
  for (j in 1:i){
    C[i,j] = sqrt(sum((CMD[i,]-CMD[j,])^2))
    D[i,j] = (C[i,j]-Y[i,j])^2
  }
}

#Average Position
PosMean = array(0, c(n,d))
for (i in 1:n){
  for (t in 1:d){
    PosMean[i,t] = mean(S[i,t,])
  }
}

#Plots
points(True,col="black",pch = 20)
points(S[,,which.min(M)],col="black",pch = 4)
points(PosMean, pch = 9)
CMDstar = procrustes(True, CMD, scale = F)
points(CMDstar$Yrot + cbind(rep(CMDstar$translation[,1],n),rep(CMDstar$translation[,2],n)), col = "red", pch = 9)
text(CMDstar$Yrot + cbind(rep(CMDstar$translation[,1],n),rep(CMDstar$translation[,2],n)), labels = 1:n, col="red")

#Compare to "true" positions
# d bar
E = array(0, c(n,n))
F = array(0, c(n,n))
G = array(0, c(n,n))
for (i in 1:n){
  for (j in 1:i){
    E[i,j] = sqrt(sum((True[i,]-True[j,])^2))
    F[i,j] = (C[i,j]-E[i,j])^2 #True vs CMD
    G[i,j] = (E[i,j]-Y[i,j])^2 #True vs average
    
  }
}

#Compare scores
min(ML) #"best" vs all other samples
sum(B) #CMD vs all samples
min(M) #"best" vs average distance Y #
sum(D) #CMD vs average distance #
sum(F) #True vs CMD
sum(G) #True vs average distance
#PosMean vs average distance #
