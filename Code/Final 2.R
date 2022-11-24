#################################################
rm(list=ls())
load("~/OneDrive/UCD/Research/GaussianLatentSpace.RData")
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
score1 = min(M) #"best" vs Y

end_time <- Sys.time()
time1 = end_time - start_time

#################################################
rm(list=ls())
load("~/Desktop/GaussianLatentSpace.RData")
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
score2 = sum(D) #CMD vs Y

end_time <- Sys.time()
time2 = end_time - start_time

#################################################
rm(list=ls())
load("~/Desktop/GaussianLatentSpace.RData")
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
score3 = sum(G) #PosMean vs Y

end_time <- Sys.time()
time3 = end_time - start_time

#################################################
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
