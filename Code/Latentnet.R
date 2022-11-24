library(latentnet)

data(samplk)
ls()
samplk.tot.m<-as.matrix(samplk1)+as.matrix(samplk2)+as.matrix(samplk3)
samplk.tot <- as.network(samplk.tot.m, directed=TRUE, matrix.type="a",
                         ignore.eval=FALSE, names.eval="nominations" # Important!
                         )
samplk.tot %v% "group" <- samplk1 %v% "group"

as.matrix(samplk.tot,attrname="nominations")[1:5,1:5] #Sociomatrix
as.matrix(samplk.tot)[1:5,1:5] #Edges
samplk.ecol <- gray(1 - (samplk.tot %e% "nominations")/4)
plot(samplk.tot,displaylabels=TRUE,
     vertex.col=as.numeric(factor(samplk.tot %v% "group"))+1,
     edge.col=samplk.ecol)

#Fixed Effects
samplk.nm.l <- ergmm(samplk.tot~nodematch("group",diff=TRUE),tofit="mle", verbose=TRUE)
samplk.nm.e <- ergm(samplk.tot~edges+nodematch("group",diff=TRUE))
summary(samplk.nm.l)
summary(samplk.nm.e)

#Latent Position Effects
samplk.d2<-ergmm(samplk.tot~euclidean(d=2), verbose=TRUE,tofit="mle")
plot(samplk.d2,what="mle", pie=TRUE,labels=TRUE)
samplk.d2G3<-ergmm(samplk.tot~euclidean(d=2,G=3), verbose=TRUE)
mcmc.diagnostics(samplk.d2G3)
plot(samplk.d2G3, pie=TRUE,labels=TRUE)

#Random Effects
samplk.d2G3r<-ergmm(samplk.tot~euclidean(d=2,G=3)+rreceiver, verbose=TRUE)
mcmc.diagnostics(samplk.d2G3r)
par(mfrow=c(1,2))
Z.K.ref <- summary(samplk.d2G3,point.est="pmean")$pmean$Z.K
Z.ref <- plot(samplk.d2G3, pie=TRUE, Z.K.ref=Z.K.ref)
plot(samplk.d2G3r, rand.eff="receiver", pie=TRUE, Z.ref=Z.ref, Z.K.ref=Z.K.ref)

#Ergmm Fits
names(samplk.d2G3r)
names(summary(samplk.d2G3r))
names(summary(samplk.d2G3r)$pmean)

#Prediction & Simulation
(pred <- predict(samplk.d2G3r))[1:5,1:5]
heatmap(predict(samplk.d2G3r),Rowv=NA,Colv=NA)
par(mfrow=c(1,3))
for(i in 1:3){
  plot(simulate(samplk.d2G3r),displaylabels=TRUE)}

#Optimal number of clusters
fits<-lapply(1:4, function(G){
  ergmm(samplk.tot~euclidean(d=2,G=G))
  }
)
do.call(rbind,lapply(lapply(fits,bic.ergmm),as.data.frame))
plot(1:4,sapply(lapply(fits,bic.ergmm),"[[","overall"),xlab="G",ylab="BIC",main="")

#Plots
plot(samplk.d2G3r, what="cloud", rand.eff="receiver", Z.ref=Z.ref, Z.K.ref=Z.K.ref)
plot(samplk.d2G3r, what="density", rand.eff="receiver", Z.ref=Z.ref, Z.K.ref=Z.K.ref)
par(mfrow=c(1,1))
for(i in 1:samplk.d2G3r$control$sample.size){
  plot(samplk.d2G3r, what=i, rand.eff="receiver")
  Sys.sleep(0.1)
}
samplk.d3G3r <- ergmm(samplk.tot~euclidean(d=3,G=3)+rreceiver)
plot(samplk.d3G3r, rand.eff="receiver", use.rgl=TRUE, labels=TRUE, main="", xlab="", ylab="")
