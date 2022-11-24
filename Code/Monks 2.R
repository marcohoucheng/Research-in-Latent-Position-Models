#################################################
rm(list=ls())
load("~/OneDrive/UCD/Research/Monks.RData")

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

#Plot the positions
par(mfrow=c(1,1))
mklpro = procrustes(OS, mkl, scale = F)
mklstar = mklpro$Yrot + cbind(rep(mklpro$translation[,1],n),rep(mklpro$translation[,2],n))
pmodepro = procrustes(OS, pmode, scale = F)
pmodestar = pmodepro$Yrot + cbind(rep(pmodepro$translation[,1],n),rep(pmodepro$translation[,2],n))
pmeanpro = procrustes(OS, pmean, scale = F)
pmeanstar = pmeanpro$Yrot + cbind(rep(pmeanpro$translation[,1],n),rep(pmeanpro$translation[,2],n))
colour = distinctColorPalette(n)
plot(CMDstar, col = colour, pch = 17, cex = 3, main = "CMD", xlab = "", ylab = "", cex.main = 1.5, xlim = c(-3,3), ylim = c(-3,4), xaxt = 'n', yaxt = 'n')
text(CMDstar, col = colour, labels = Names, pos = 4, cex = 2)
grid()
plot(pmodestar, col = colour, pch = 19, cex = 3, main = "Posterior Mode", xlab = "", ylab = "", cex.main = 1.5, xlim = c(-3,3), ylim = c(-3,4), xaxt = 'n', yaxt = 'n')
text(pmodestar, col = colour, labels = Names, pos = 4, cex = 2)
grid()
plot(pmodestar, col = colour, pch = 19, cex = 2.5, xlab = "", ylab = "", xlim = c(-3,3), ylim = c(-3,4), xaxt = 'n', yaxt = 'n', ann = FALSE)
#text(pmodestar, col = colour, labels = Names, pos = 4, cex = 1.5)
arrows(x0 = pmodestar[,1], y0 = pmodestar[,2], x1 = CMDstar[,1], y1 = CMDstar[,2], col = colour, lwd = 3)
grid()

legend("topright", inset = 0.02, legend = c("CMD", "Posterior Mode"), pch = c(19, 17))

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

scorebar <- barplot(average_edges, names.arg = c("Bayes Estimator", "CMD", "Average Position", "MKL", "Posterior Mode"), main = "Expected Number of Total Edges", ylim = c(0,100))
text(scorebar,average_edges+2,labels=as.character(round(average_edges)))
abline(h=True_no_of_edges, col = "Red", lty = 5, lwd = 2.5)
