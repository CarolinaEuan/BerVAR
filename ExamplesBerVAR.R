### Examples of BerVAr
#Example 1 

source('SourceBerVAR.R')
m<-3
states<-as.binary((1:2^m)-1,n=m,littleEndian=TRUE)
states.names<-(do.call(rbind, states))*1
states.names<-paste0(states.names[,1],states.names[,2],states.names[,3])
states.names

#Marginals pi and pi^C
m<-3
M<-2^m
##Same marginals
p_y<-c(.7,.7,.7)
p_z<-c(.4,.3,.4)

Omega<-matrix(as.numeric(unlist(as.binary((1:M)-1,n=m,littleEndian=TRUE))),ncol=m,byrow = TRUE)

pho_y<-c(-.3,-.3,.3)
pho_z<-c(-.6,.2,.2)
#Different covariance structure
#Multivariate bernoulli densities
pm_y<-Generate.pi.Mod(pho_y,p_y,MOD=0,Omega)
pm_z<-Generate.pi.Mod(pho_z,p_z,MOD=0,Omega)


par(mfrow=c(2,1))
barplot(t(pm_y),names.arg = states.names,legend.text = FALSE,ylim=c(0,.5),
        main=expression(paste("Distribution of ",Y[t])),xlab="States",ylab="Probability")
barplot(t(pm_z),names.arg = states.names,legend.text = FALSE,ylim=c(0,.5),
        main=expression(paste("Distribution of ",Z[t])),
        xlab="States",ylab="Probability")



#Time series
#3914
set.seed(3914)
X0<-rbinom(m,1,.5)
TT<-2000+50
Xt<-matrix(NA,ncol=m,nrow=TT)
#Simulate latent states 
ksample_y<-sample(1:(2^m),prob = pm_y,replace = TRUE,size = TT)
ksample_z<-sample(1:(2^m),prob = pm_z,replace = TRUE,size = TT)

Xt[1,]<-Omega[ksample_y[1],]* X0 + Omega[ksample_z[1],]*(1 - X0)

for(i in 1:(TT-1)) {
  Xt[i+1,]<-Omega[ksample_y[1+i],]* Xt[i,] + 
    Omega[ksample_z[1+i],]*(1 - Xt[i,]) }

Xt<-Xt[51:TT,]


#Xt is the simulated TS
layout(matrix(c(1,1,1,2,2,2,3,3,3,3),ncol=1))
par(mar=c(0,5,1,1))
barplot(Xt[200:300,1],axes=FALSE,ylab = expression(X[paste(1,t)]),cex.lab=2,col=4)
axis(2,at=c(0,1),c(0,1),cex.axis=1.5)
barplot(Xt[200:300,2],axes=FALSE,ylab = expression(X[paste(2,t)]),cex.lab=2,col=2)
axis(2,at=c(0,1),c(0,1),cex.axis=1.5)
par(mar=c(4,5,1,1))
barplot(Xt[200:300,3],axes=FALSE,ylab = expression(X[paste(3,t)]),cex.lab=2,xlab = "Time")
axis(2,at=c(0,1),c(0,1),cex.axis=1.5)
axis(1,cex.axis=1.5)

mu<-p_z/(1-p_y+p_z)
sigma<-((1-p_y)*p_y*(mu^2)+(1-p_z)*p_z*(1-mu)^2)/(p_y+p_z-2*p_y*p_z)

#Mean and Variance Sample vs True
rbind(mu,apply(Xt, 2, "mean"),sigma,apply(Xt, 2, "var"))
round(rbind(mu,apply(Xt, 2, "mean"),sigma,apply(Xt, 2, "var")),3)

#Autocovariance and Crosscovariance - Sample Vs True
Gamma<-acf(Xt,plot=FALSE,lag.max = 15,type = "correlation")
#CovEmp<-function(h,i,j)cor(Xt[1:(2000-h),i], Xt[h:2000,j], method = "pearson")

GammaMat<-CovMat(p_y,p_z,pho_y,pho_z,max.lag=15)
indexij<-cbind(rep(1:3,3),rep(1:3,each=3))
normconst<-GammaMat[1,]
normconst[2:3]<-sqrt(GammaMat[1,1])*sqrt(c(GammaMat[1,5],GammaMat[1,9]))
normconst[c(4,6)]<-sqrt(GammaMat[1,5])*sqrt(c(GammaMat[1,1],GammaMat[1,9]))
normconst[c(7,8)]<-sqrt(GammaMat[1,9])*sqrt(c(GammaMat[1,1],GammaMat[1,5]))

par(mfrow=c(3,3))
for(i in 1:9){
  plot(0:15,Gamma$acf[,indexij[i,1],indexij[i,2]],ylim=c(min(Gamma$acf)-.1,1),
       ylab = " ",xlab="lag",pch=20,col="gray",
       main = mymain(indexij[i,1],indexij[i,2]))
  lines(0:15,GammaMat[,i]/normconst[i],col="blue")
  abline(h=0,col=2,lty=2)
}


#MLE Estimation 
#Full Model
Est.QMLE<-try(MLE.BerVAR(t(Xt),Omega,MOD = 0),silent = TRUE)
#Parsimonious
Est.QMLE.S<-try(MLE.BerVAR(t(Xt),Omega,MOD = 1),silent=TRUE)

list(Est.QMLE,Est.QMLE.S)

#Observed estimated laws
#Full MOD
Fpiy.hat<-Generate.pi.Mod(par = Est.QMLE[[2]]$Full.model.par.y[,1],p = Est.QMLE[[2]]$p_y[,1],MOD=0,Omega)
Fpiz.hat<-Generate.pi.Mod(par = Est.QMLE[[4]]$Full.model.par.z[,1],p = Est.QMLE[[4]]$p_z[,1],MOD=0,Omega)

#Parsimonnious
Ppiy.hat<-Generate.pi.Mod(par = Est.QMLE.S[[2]]$Single.model.par.y[,1],p = Est.QMLE.S[[2]]$p_y[,1],MOD=1,Omega)
Ppiz.hat<-Generate.pi.Mod(par = Est.QMLE.S[[4]]$Single.model.par.z[,1],p = Est.QMLE.S[[4]]$p_z[,1],MOD=1,Omega)


par(mfrow=c(1,2))
plot(pm_y,pch=20,main = expression(paste("Estimated Density of ",Y[t])), xlab=expression(Omega[M]),
     ylab ="Probability",axes=FALSE ,ylim = c(0,max(c(pm_z,Fpiy.hat))+.2))
points(Fpiy.hat,pch=2,col=4)
points(Ppiy.hat,pch=3,col=2)
axis(2)
box()
axis(1,at = 1:M,labels = states.names)
legend("topright",legend = c("True","Full Model","Parsimonious Model"),pch=c(20,2,3),col=c(1,4,2),bty="n")


plot(pm_z,pch=20,main = expression(paste("Estimated Density of ",Z[t])), xlab=expression(Omega[M]),
     ylab ="Probability",axes=FALSE ,ylim = c(0,max(c(pm_z,Ppiz.hat))+.2))
points(Fpiz.hat,pch=2,col=4)
points(Ppiz.hat,pch=3,col=2)
axis(2)
box()
axis(1,at = 1:M,labels = states.names)
legend("topright",legend = c("True","Full Model","Parsimonious Model"),pch=c(20,2,3),col=c(1,4,2),bty="n")



TCov<-CovMat(pY=p_y,pZ=p_z,phoY = pho_y,phoZ = pho_z,max.lag = 10)
FCov<-CovMat(pY=Est.QMLE[[2]]$p_y[,1],pZ=Est.QMLE[[4]]$p_z[,1],
             phoY =Est.QMLE[[2]]$Full.model.par.y[,1],phoZ = Est.QMLE[[4]]$Full.model.par.z[,1],max.lag = 10)
PCov<-CovMat(pY=Est.QMLE.S[[2]]$p_y[,1],pZ=Est.QMLE.S[[4]]$p_z[,1],
             phoY =Est.QMLE.S[[2]]$Single.model.par.y[,1],phoZ = Est.QMLE.S[[4]]$Single.model.par.z[,1],max.lag = 10)

indexij<-cbind(rep(1:m,m),rep(1:m,each=m))
ind.equal<-which(indexij[,1]==indexij[,2])
gamma0<-TCov[1,ind.equal]
Fgamma0<-FCov[1,ind.equal]
Pgamma0<-PCov[1,ind.equal]

#Estimated correlation plots
  par(mfrow=c(m,m))
  for(i in 1:(m^2)){
    plot(0:10,FCov[,i]/sqrt(prod(Fgamma0[indexij[i,]])),ylim=c(-.2,1),
         ylab = " ",xlab="lag",col="blue",type="l",lwd=2,lty=2,
         main = mymain(indexij[i,1],indexij[i,2]))
    lines(0:10,PCov[,i]/sqrt(prod(Pgamma0[indexij[i,]])),col=2,lty=3,lwd=2)
    lines(0:10,TCov[,i]/sqrt(prod(gamma0[indexij[i,]])))
    abline(h=0,col="gray",lty=2)
  }




