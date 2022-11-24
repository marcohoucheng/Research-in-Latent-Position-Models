rm(list=ls())

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
  plot(S[,,1], col=colour, xlab="", ylab="")
  for (i in 2:k){
    points(S[,,i],col=colour)
  }
  assign("S", S, envir=globalenv())
  assign("True", True, envir=globalenv())
}

GaussianLatentSpaceGenerator(n,d,mu,var,k,Var)
# T true data
# S simulated data


##############################################################################################################
#Original

M = array(0, dim=k)
I = array(0, dim=c(k,n))
A = array(0, dim=c(k,n,n)) #This is the extra term that we don't need
J = array(0, dim=c(k,n,n))
X = array(0, dim=c(k,n,n))
z = array(0, dim=c(k,n,n))
K = array(0, dim=c(k,n,n,k))
Y = array(0, dim=c(k,n,n,k))
x = array(0, dim=c(k,n,n,d))
y = array(0, dim=c(k,n,n,k,d))
for (m in 1:k){
  #I = array(0, dim=c(m,n))
  for (i in 1:n){
    #J = array(0, dim=c(m,n,n))
    #X = array(0, dim=c(m,n,n))
    #z = array(0, dim=c(m,n,n))
    for (j in 1:n){ #Change n to i
      #K = array(0, dim=c(m,i,j,k))
      #Y = array(0, dim=c(m,i,j,k))
      #x = array(0, dim=c(m,i,j,d))
      for (t in 1:k){
        #y = array(0, dim=c(m,i,j,t,d))
        for (f in 1:d){
          y[m,i,j,t,f] = (S[i,f,t]-S[j,f,t])^2
        }
        Y[m,i,j,t] = sum(y[m,i,j,t,])
      }
      A[m,i,j] = sum(Y[m,i,j,])
      z[m,i,j] = sum(sqrt(Y[m,i,j,]))
      #for (f in 1:d){
        #x[m,i,j,f] = (S[i,f,m]-S[j,f,m])^2
      #}
      #x[m,i,j,] = y[m,i,j,m,]
      #X[m,i,j] = sum(x[m,i,j,])
      #X[m,i,j] = sum(y[m,i,j,m,])
      X[m,i,j] = Y[m,i,j,m]
      J[m,i,j] = k*X[m,i,j]-2*sqrt(X[m,i,j])*z[m,i,j]+A[m,i,j]
    }
    I[m,i] = sum(J[m,i,])
  }
  M[m] = sum(I[m,])
}
min(M)
which.min(M)
M[which.min(M)]
points(S[,,which.min(M)],col="black",pch=2)

#How big is A?
A <- A[1,,]
A

##############################################################################################################
#1
start_time <- Sys.time()
M = array(0, k)
X = array(0, c(n,n))
Z = array(0, c(n,n,k))
for (m in 1:k){
  for (i in 1:n){
    for (j in 1:i){
      for (t in 1:k){
        Z[i,j,t] = sum((S[i,,t]-S[j,,t])^2)
      }
      X[i,j] = k*Z[i,j,m]-2*sqrt(Z[i,j,m])*sum(sqrt(Z[i,j,]))
    }
  }
  M[m] = sum(X)
}
min(M)
which.min(M)
end_time <- Sys.time()
end_time - start_time

M[which.min(M)]
points(S[,,which.min(M)],col="black",pch=2)

##############################################################################################################
#2
start_time <- Sys.time()
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

min(ML)
which.min(ML)
end_time <- Sys.time()
end_time - start_time

M[which.min(ML)]
points(S[,,which.min(ML)],col="black",pch=4)

##############################################################################################################
#3
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
    Y[i,j] = k^(-1)*sum(Z[i,j,]) #d bar
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

min(M)
which.min(M)
end_time <- Sys.time()
end_time - start_time

M[which.min(M)]
points(S[,,which.min(M)],col="black",pch=4)

##############################################################################################################
D = array(0, c(n,n))
D = Y + t(Y)
CMD = cmdscale(D, k = 2)
points(CMD,col="black",pch=15)

SAM = sammon(D, k = 2)
points(SAM$points,col="black",pch=19)

# Can fix the position problem, but how do we know which node represents which actor?
CMDbar = c(mean(CMD[,1]), mean(CMD[,2]))
drift = c(mean(S[,1,which.min(M)]), mean(S[,2,which.min(M)]))
Z0 = S[,,which.min(M)] - drift #Z0 centred at 0
Z1 = CMD - CMDbar #transform Z centred at 0
Zstar = Z0%*%t(Z1)%*%Re(solve(sqrtm(Z1%*%t(Z0)%*%Z0%*%t(Z1))))%*%Z1
CMDstar = Zstar + drift
points(CMDstar, col = "black", pch = 9)

CMDstar = procrustes(S[,,which.min(M)], CMD, scale = FALSE, symmetric = FALSE)$Yrot
points(CMDstar, col = "black", pch = 9)


# Compare S[,,which.min(M)] (done) to CMDstar in terms of the 2 loss functions

# optimal sample

A = array(0, c(n,n))
B = array(0, c(n,n))

for (i in 1:n){
  for (j in 1:i){
    A[i,j] = sum((CMDstar[i,]-CMDstar[j,])^2)
    B[i,j] = k*A[i,j]-2*sqrt(A[i,j])*YL[i,j]
  }
}
sum(B)

# d bar
C = array(0, c(n,n))
D = array(0, c(n,n))
for (i in 1:n){
  for (j in 1:i){
    C[i,j] = sqrt(sum((CMDstar[i,]-CMDstar[j,])^2))
    D[i,j] = (C[i,j]-Y[i,j])^2
  }
}
sum(D)

# the optimal sample from the iterations is better than the CMD...


## Jordan
Re(solve(sqrtm(Z1%*%t(Z0)%*%Z0%*%t(Z1))))

TI = solve(Z1%*%t(Z0)%*%Z0%*%t(Z1))
E = eigen(TI)
V = E$vectors
Jordan = Diagonal( x = E$values )

Vinv = solve(V)
V %*% J %*% Vinv # = TI

Jsqrt = Diagonal(x=sqrt(round(E$values)))
Jsqrt = Diagonal(x=sqrt(round(E$values, 5)))
V %*% Jsqrt %*% Vinv # = TI^0.5
inv = V %*% Jsqrt %*% Vinv
inv = solve(V %*% Jsqrt %*% Vinv)

