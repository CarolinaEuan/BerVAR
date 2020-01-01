#############################################
# Source File to fit the BerVAR model to data
#############################################

#Required libraries for numeric computations
library(fields)
library(fda)
library(car)
library(binaryLogic)
library(lpSolve)
library(bigstatsr)
library(Matrix)
library(klin)
library(limSolve)
library(numDeriv)
#require sudo apt-get install libgmp3
library(rcdd)

#Defined functions: 

# mymain - auxiliarity function to create plot title in CovMat
# CovMat - computes the covariance function of the BerVAR process guven the parameters of the latent process
# Generate.pi - generates pi of the multivariate Bernoulli density by considering Teuges formula 
# Check.pho - finds valid values for the correlation parameters using the inequalities shown in Fontana and Semeraro 
# Aux.R - given p finds the ray densities for pi (see Fontana and Semeraro)
# Find.joint - given p and correlation parameters finds a valid multivariate Bernoulli 
#              density by following Fontana and Semarano procedure
# MLE.BerVAR - fits the BerVAR model to a multivariate binary TS X using the proposed method
# Generate.pi.Mod - given p and correlation parameters and assuming one of the proposed submodels 
#                   (Full, parsimonious or Spatial) computes the valid bernoulli density
# PredPi - Computes the predictive distribution of X_{t+h} | X_{t}


mymain<-function(k,l,Corr=TRUE){
  if(Corr){
  if(k==l){return(bquote(paste("ACF of ",X[paste(.(k),t)])))
  }else{return(bquote(paste("Correlation of ",X[paste(.(k),t)]," and ",X[paste(.(l),t)])))}
  }else{
    if(k==l){return(bquote(paste("Covariance of ",X[paste(.(k),t)])))
    }else{return(bquote(paste("Covariance of ",X[paste(.(k),t)]," and ",X[paste(.(l),t)])))}
  }
}

CovMat<-function(pY,pZ,phoY,phoZ,max.lag=15){
  #Covariance matrix by col
  #row is the lag and col is the index
  #index sequence is 1,1 1,2 1,3....1,m, 2,1, 2,2...
  h<-seq(0,max.lag,by = 1)
  mu<-pZ/(1-pY+pZ)
   phoYMat<-matrix(1,ncol=length(pY),nrow =length(pY))
   phoYMat[upper.tri(phoYMat)]<-phoY
   phoYMat[lower.tri(phoYMat)]<-phoY
   phoZMat<-matrix(1,ncol=length(pZ),nrow =length(pZ))
   phoZMat[upper.tri(phoZMat)]<-phoZ
   phoZMat[lower.tri(phoZMat)]<-phoZ
  Gam<-matrix(0,nrow=length(h),ncol=length(pY)^2)
  cont<-0
  for(i in 1:length(pY))for(j in 1:length(pY)){
    cont<-cont+1
    if(i==j){
      gamii0<-((1-pY[i])*pY[i]*mu[i]^2+(1-pZ[i])*pZ[i]*(1-mu[i])^2)/(1-(pY[i]+pZ[i]-2*pY[i]*pZ[i]))
      gam<-((pY[i]-pZ[i])^h)*gamii0 
    }else{
      thetaYij<-phoYMat[i,j]*sqrt(pY[i]*(1-pY[i]))*sqrt(pY[j]*(1-pY[j]))
      thetaZij<-phoZMat[i,j]*sqrt(pZ[i]*(1-pZ[i]))*sqrt(pZ[j]*(1-pZ[j]))
      gamij0<-(thetaZij*(1-mu[i]-mu[j])+mu[i]*mu[j]*(thetaYij+thetaZij))/
              (1-(thetaYij+thetaZij)-(pY[i]-pZ[i])*(pY[j]-pZ[j]))
      gam<-((pY[j]-pZ[j])^h)*gamij0
    }
    Gam[,cont]<-gam
  }
  colnames(Gam)<-paste("Cov",rep(1:length(pY),each=length(pZ)),rep(1:length(pZ),length(pY)))
  return(Gam)
}

Generate.pi<-function(p,s)
{
  #m>1
  #Generates pi by considering Teuges formula
  #all moments s need to be given
  m<-length(p)
  A <- lapply(rev(1:m),function(i) Matrix(c(1-p[i],-1,p[i],1),2,2,byrow= TRUE))
  pi<-klin.eval(A, s) 
  return(pi)
}

Check.pho<-function(p,Omega){
  m<-length(p)
  M<-2^m
  #Omega<-t(sapply(0:(M-1),function(x){ as.integer(intToBits(x))}))[,1:m]
  SD<-sqrt(p*(1-p))
  Order.Moment<-rowSums(Omega)
  #Compute R matrix
  Order0<-which(Order.Moment==1)
  Rp<-Aux.R(p,Order0)
  #Compute bounds
  Order2<-which(Order.Moment==2)
  B <- lapply(rev(1:m),function(i) t(Matrix(c(1,0,1,1),2,2)))
  A2<-matrix(0,ncol=M,nrow = length(Order2))
  for(i in 1:length(Order2)){
    ei<-rep(0,M)
    ei[Order2[i]]<-1
    Ai<-klin.eval(B,ei)
    A2[i,]<-Ai
  }
  Bounds<-matrix(0,ncol=2,nrow = length(Order2))
  A2alpha<-A2%*%Rp
  for(i in 1:length(Order2)){
    Bounds[i,]<-c((min(A2alpha[i,])-prod(p[Omega[Order2[i],]==1]))/prod(SD[Omega[Order2[i],]==1]),
                  (max(A2alpha[i,])-prod(p[Omega[Order2[i],]==1]))/prod(SD[Omega[Order2[i],]==1]))
  }
  return(Bounds)
}

Aux.R<-function(p,Order0){
  m<-length(p)
  M<-2^m
  A <- lapply(rev(1:m),function(i) t(Matrix(c(1,1,-p[i],1-p[i]),2,2,byrow= TRUE)))
  AO<-matrix(0,ncol=M,nrow = length(Order0))
  for(i in 1:length(Order0)){
    ei<-rep(0,2^m)
    ei[Order0[i]]<-1
    Ai<-klin.eval(A,ei)
    AO[i,]<-Ai
  }
  if(m<6){
  II<-diag(1,M,M)
  MM<-makeH(a1 = -II,b1 = rep(0,M),a2 =rbind(AO,rep(1,M)),b2 = c(rep(0,nrow(AO)),1) )
  VV<-scdd(MM,representation = "H")
  Rp<-t(VV$output[,-c(1:2)])
  }else{
    Samples<-xsample(E=rbind(AO,matrix(1,ncol=ncol(AO),nrow=1)),F=c(rep(0,m),1),
                     G=diag(ncol(AO)),H = rep(0,ncol(AO)),outputlength = m*(m+1)+10)
    Rp<-t(Samples$X)
  }
  if(ncol(Rp)>200) Rp<-round(Rp[,1:200],10)
  return(Rp)
}

Find.joint<-function(p,pho,Omega,eps=exp(-15)){
  m<-length(p)
  M<-2^m
  Order.Moment<-rowSums(Omega)
  #Second order moments mu_2
  Order2<-which(Order.Moment==2)
  SD<-sqrt(p*(1-p))
  M2<-pho
  for(i in 1:length(pho)){
    #if(pho[i]==0)M2[i]<-prod(p[Omega[Order2[i],]==1])else 
    M2[i]<-pho[i]*prod(SD[Omega[Order2[i],]==1])+prod(p[Omega[Order2[i],]==1])
  }
  #Compute M matrix M%*%pi=mu
  B <- lapply(rev(1:m),function(i) t(Matrix(c(1,0,1,1),2,2)))
  #Submatrix of order 2
  A2<-matrix(0,ncol=M,nrow = length(Order2))
  for(i in 1:length(Order2)){
    ei<-rep(0,M)
    ei[Order2[i]]<-1
    Ai<-klin.eval(B,ei)
    A2[i,]<-Ai
  }
  #If m is small we compute the ray matrix
  #Compute R matrix
  Order0<-which(Order.Moment==1)
  Rp<-Aux.R(p,Order0)
  A2alpha<-A2%*%Rp  
  lambdaAux<-lsei(A=A2alpha,B=M2,E=matrix(1,ncol=ncol(A2alpha),nrow=1),F=1,
                  G=diag(ncol(A2alpha)),H = rep(eps,ncol(A2alpha)))
  lambda<-lambdaAux$X
  pi<-Rp%*%lambda
  return(pi)
}

MLE.BerVAR<-function(X,Omega,MOD=0,loc=NULL,bige=exp(15),smalle=exp(-15),par0z=NULL,par0y=NULL){
  #X = observed time series by row
  m<-nrow(X)
  M<-2^m
  TT<-ncol(X)
  Order.Moment<-rowSums(Omega)
  Order0<-which(Order.Moment==1)
  Xpast<-X[,1:(TT-1)]
  Xt<-X[,2:TT]
  FIT<-vector("list",4)
  if(MOD==2) dist<-apply(Omega[Order.Moment==2,], 1, f<-function(x) sqrt(diff(loc[x==1,1])^2+diff(loc[x==1,2])^2))
  #Likelihood maximizes marginally
  #Y
  FIT[[1]]<-"Y Process fitting QMLE results"
  mod.names<-c("Full","Single","Spatial")
  FIT[[2]]<-vector("list",2)
  names(FIT[[2]])<-c("p_y",paste0(mod.names[MOD+1],".model.par.y"))
  phat_y<-rowSums(Xpast*Xt)/rowSums(Xpast==1)
  ny<-rowSums(Xpast==1)
  Sy<-rowSums(Xpast*Xt)
  var.mle<-Sy*(ny-Sy)/ny^3
  Results.p<-cbind(phat_y,phat_y-2*sqrt(var.mle),phat_y+2*sqrt(var.mle))
  colnames(Results.p)<-c("MLE","LI","UI")
  FIT[[2]]$p_y<-Results.p
  #Compute R matrix
  Ry<-Aux.R(phat_y,Order0)
  SDy<-sqrt(phat_y*(1-phat_y))
  
  Order2<-which(Order.Moment==2)
  B <- lapply(rev(1:m),function(i) t(Matrix(c(1,0,1,1),2,2)))
  A2<-matrix(0,ncol=M,nrow = length(Order2))
  for(i in 1:length(Order2)){
    ei<-rep(0,M)
    ei[Order2[i]]<-1
    Ai<-klin.eval(B,ei)
    A2[i,]<-Ai
  }

  Pse.Log.Lik.Y<-function(par){
    #We choose a model  
    if(MOD==0)pho<-par
    if(MOD==1)pho<-rep(par,m*(m-1)/2)
    if(MOD==2){
      #If MOD = 2 par=lambda) and loc != NULL
      lam<-par #a<-par[1]
      pho<-exp(-lam*dist)
      #loglik
      A2alpha<-A2%*%Ry  
      M2<-pho
      for(i in 1:length(pho)){
        M2[i]<-pho[i]*prod(SDy[Omega[Order2[i],]==1])+prod(phat_y[Omega[Order2[i],]==1])
      }
      lambdaAux<-lsei(A=A2alpha,B=M2,E=matrix(1,ncol=ncol(A2alpha),nrow=1),F=1,
                      G=diag(ncol(A2alpha)),H = rep(smalle,ncol(A2alpha)))
      lambda<-lambdaAux$X
      pi<-Ry%*%lambda
      logl<-0
      if(any(pi < 0)){logl<- -bige}else{
        log1<-0
        for(i in 1:(TT-1)){
          y.aux<-Xt[Xpast[,i]==1,i]
          if(length(y.aux)==0){ logl1<-0 }else{
            if(length(y.aux)==1){
              logl1<-log(sum(pi[Omega[,Xpast[,i]==1]==y.aux]))
            }else{
              row.is.a.match <- apply(Omega[,Xpast[,i]==1], 1, all.equal, current=y.aux) 
              match.idx <- which(row.is.a.match==TRUE) 
              logl1<-log(sum(pi[match.idx]))  
            }
          }
          logl<-ifelse(logl1==-Inf,logl,logl+logl1)
          #print(c(logl,i))
        }
      }
      return(logl)
    }
    #We check model if feasible (if MOD!=2)
    Bounds<-matrix(0,ncol=2,nrow = length(Order2))
    A2alpha<-A2%*%Ry
    for(i in 1:length(Order2)){
      Bounds[i,]<-c((min(A2alpha[i,])-prod(phat_y[Omega[Order2[i],]==1]))/prod(SDy[Omega[Order2[i],]==1]),
                    (max(A2alpha[i,])-prod(phat_y[Omega[Order2[i],]==1]))/prod(SDy[Omega[Order2[i],]==1]))
    }
    if((sum(pho<Bounds[,1],na.rm = TRUE)+sum(pho>Bounds[,2],na.rm = TRUE))>0){
      return(-bige)
    }else{
      A2alpha<-A2%*%Ry  
      M2<-pho
      for(i in 1:length(pho)){
        M2[i]<-pho[i]*prod(SDy[Omega[Order2[i],]==1])+prod(phat_y[Omega[Order2[i],]==1])
      }
      lambdaAux<-lsei(A=A2alpha,B=M2,E=matrix(1,ncol=ncol(A2alpha),nrow=1),F=1,
                      G=diag(ncol(A2alpha)),H = rep(0,ncol(A2alpha)))
      lambda<-lambdaAux$X
      pi<-Ry%*%lambda
      logl<-0
      if(any(pi < 0)){logl<- -bige}else{
        log1<-0
        for(i in 1:(TT-1)){
          y.aux<-Xt[Xpast[,i]==1,i]
          if(length(y.aux)==0){ logl1<-0 }else{
            if(length(y.aux)==1){
              logl1<-log(sum(pi[Omega[,Xpast[,i]==1]==y.aux]))
            }else{
              row.is.a.match <- apply(Omega[,Xpast[,i]==1], 1, all.equal, current=y.aux) 
              match.idx <- which(row.is.a.match==TRUE) 
              logl1<-log(sum(pi[match.idx]))  
            }
          }
          logl<-logl+logl1
        }
      }
      return(logl)  
    }
  }
  if(MOD==1){
    OPT.Y<-optimize(Pse.Log.Lik.Y, c(-1, 1), tol = 0.0001,maximum = TRUE)
    var.mle<- -1/(numDeriv::hessian(func=Pse.Log.Lik.Y, x=OPT.Y$maximum))
    Results.par<-cbind(OPT.Y$maximum,OPT.Y$maximum-2*sqrt(var.mle),OPT.Y$maximum+2*sqrt(var.mle))
  }else{
  if(MOD==0){
    if(is.null(par0y))par0<-cor(-1*(t(X)-1))[upper.tri(cor(t(X)))] else par0<-par0y
    OPT.Y<-optim(fn=Pse.Log.Lik.Y,par = par0,method = "Nelder-Mead",control=list(fnscale=-1),hessian = TRUE)
    var.mle<- -1/diag(OPT.Y$hessian)
    Results.par<-cbind(OPT.Y$par,OPT.Y$par-2*sqrt(var.mle),OPT.Y$par+2*sqrt(var.mle))
  }
  if(MOD==2){
    if(is.null(par0y))par0<-c(0.01, 10) else par0<-par0y
    OPT.Y<-optimize(Pse.Log.Lik.Y, par0, tol = 0.0001,maximum = TRUE)
    var.mle<- -1/hessian(func=Pse.Log.Lik.Y, x=OPT.Y$maximum)
    Results.par<-cbind(OPT.Y$maximum,OPT.Y$maximum-2*sqrt(var.mle),OPT.Y$maximum+2*sqrt(var.mle))
  }
  }
  colnames(Results.par)<-c("MLE","LI","UI")
  FIT[[2]][[2]]<-Results.par
  #Z
  FIT[[3]]<-"Z Process fitting QMLE results"
  FIT[[4]]<-vector("list",2)
  names(FIT[[4]])<-c("p_z",paste0(mod.names[MOD+1],".model.par.z"))
  pzhat<-rowSums((1-Xpast)*Xt)/rowSums(Xpast==0)
  nz<-rowSums(Xpast==0)
  Sz<-rowSums((1-Xpast)*Xt)
  var.mle<-Sz*(nz-Sz)/nz^3
  Results.p<-cbind(pzhat,pzhat-2*sqrt(var.mle),pzhat+2*sqrt(var.mle))
  colnames(Results.p)<-c("MLE","LI","UI")
  FIT[[4]]$p_z<-Results.p
  #Compute R matrix
  Rz<-Aux.R(pzhat,Order0)
  SDz<-sqrt(pzhat*(1-pzhat))
  Pse.Log.Lik.Z<-function(par){
    #We choose a model  
    if(MOD==0)pho<-par
    if(MOD==1)pho<-rep(par,m*(m-1)/2)
    if(MOD==2){
      #If MOD = 2 par=(a,lambda) and loc != NULL
      lam<-par
      #a<-par[1]
      pho<-exp(-lam*apply(Omega[Order.Moment==2,], 1, f<-function(x) sqrt(diff(loc[x==1,1])^2+diff(loc[x==1,2])^2)))
      #Loglik 
      A2alpha<-A2%*%Rz  
      M2<-pho
      for(i in 1:length(pho)){
        M2[i]<-pho[i]*prod(SDz[Omega[Order2[i],]==1])+prod(pzhat[Omega[Order2[i],]==1])
      }
      lambdaAux<-lsei(A=A2alpha,B=M2,E=matrix(1,ncol=ncol(A2alpha),nrow=1),F=1,
                      G=diag(ncol(A2alpha)),H = rep(smalle,ncol(A2alpha)))
      lambda<-lambdaAux$X
      piz<-Rz%*%lambda
      logl<-0
      if(any(piz < 0)){logl<- -bige}else{
        logl2<-0
        for(i in 1:(TT-1)){
          z.aux<-Xt[Xpast[,i]==0,i]
          if(length(z.aux)==0){logl2<-0}else{
            if(length(z.aux)==1){
              logl2<-log(sum(piz[Omega[,Xpast[,i]==0]==z.aux]))
            }else{
              row.is.a.match <- apply(Omega[,Xpast[,i]==0], 1, all.equal, current=z.aux) 
              match.idx <- which(row.is.a.match==TRUE) 
              logl2<-log(sum(piz[match.idx]))
            }
          }
          logl<-ifelse(logl2==-Inf,logl,logl+logl2)
          #print(logl)
        }
      }
      return(logl)  
    }
    #We check model if feasible
    Bounds<-matrix(0,ncol=2,nrow = length(Order2))
    A2alpha<-A2%*%Rz
    for(i in 1:length(Order2)){
      Bounds[i,]<-c((min(A2alpha[i,])-prod(pzhat[Omega[Order2[i],]==1]))/prod(SDz[Omega[Order2[i],]==1]),
                    (max(A2alpha[i,])-prod(pzhat[Omega[Order2[i],]==1]))/prod(SDz[Omega[Order2[i],]==1]))
    }
    if((sum(pho<Bounds[,1],na.rm = TRUE)+sum(pho>Bounds[,2],na.rm = TRUE))>0){
      logl<- -bige 
      return(logl) 
    }else{
      A2alpha<-A2%*%Rz  
      M2<-pho
      for(i in 1:length(pho)){
        M2[i]<-pho[i]*prod(SDz[Omega[Order2[i],]==1])+prod(pzhat[Omega[Order2[i],]==1])
      }
      lambdaAux<-lsei(A=A2alpha,B=M2,E=matrix(1,ncol=ncol(A2alpha),nrow=1),F=1,
                      G=diag(ncol(A2alpha)),H = rep(0,ncol(A2alpha)))
      lambda<-lambdaAux$X
      piz<-Rz%*%lambda
      logl<-0
      if(any(piz < 0)){logl<- -bige}else{
        logl2<-0
        for(i in 1:(TT-1)){
        z.aux<-Xt[Xpast[,i]==0,i]
        if(length(z.aux)==0){logl2<-0}else{
          if(length(z.aux)==1){
            logl2<-log(sum(piz[Omega[,Xpast[,i]==0]==z.aux]))
          }else{
            row.is.a.match <- apply(Omega[,Xpast[,i]==0], 1, all.equal, current=z.aux) 
            match.idx <- which(row.is.a.match==TRUE) 
            logl2<-log(sum(piz[match.idx]))
          }
        }
        logl<-logl+logl2 
        }
      }
      return(logl)  
    }
  }
  if(MOD==1){
    OPT.Z<-optimize(Pse.Log.Lik.Z, c(-1, 1), tol = 0.0001,maximum = TRUE)
    var.mle<- -1/(numDeriv::hessian(func=Pse.Log.Lik.Z, x=OPT.Z$maximum))
    Results.par<-cbind(OPT.Z$maximum,OPT.Z$maximum-2*sqrt(var.mle),OPT.Z$maximum+2*sqrt(var.mle))
  }else{
  if(MOD==0){
    if(is.null(par0z))par0<-cor(-1*(t(X)-1))[upper.tri(cor(t(X)))] else par0<-par0z
    OPT.Z<-optim(fn=Pse.Log.Lik.Z,par = par0,method = "Nelder-Mead",control=list(fnscale=-1),hessian = TRUE)
    var.mle<- -1/diag(OPT.Z$hessian)
    Results.par<-cbind(OPT.Z$par,OPT.Z$par-2*sqrt(var.mle),OPT.Z$par+2*sqrt(var.mle))
    }
  if(MOD==2){
    if(is.null(par0y))par0<-c(0.01, 10) else par0<-par0z
  OPT.Z<-optimize(Pse.Log.Lik.Z, par0, tol = 0.0001,maximum = TRUE)
  var.mle<- -1/hessian(func=Pse.Log.Lik.Z, x=OPT.Z$maximum)
  Results.par<-cbind(OPT.Z$maximum,OPT.Z$maximum-2*sqrt(var.mle),OPT.Z$maximum+2*sqrt(var.mle))
  }
  }
  colnames(Results.par)<-c("MLE","LI","UI")
  FIT[[4]][[2]]<-Results.par
  return(FIT)
}

Generate.pi.Mod<-function(par,p,MOD,Omega,loc=NULL,eps=exp(-15)){
  m<-length(p)
  M<-2^m
  Order.Moment<-rowSums(Omega)
  SD<-sqrt(p*(1-p))
  if(MOD==2){
      #If MOD = 2 par=lambda and loc != NULL
      lam<-par
      pho<-exp(-lam*apply(Omega[Order.Moment==2,], 1, f<-function(x) sqrt(diff(loc[x==1,1])^2+diff(loc[x==1,2])^2)))
      pi<-Find.joint(p,pho,Omega = Omega,eps=eps)
      return(pi)
  }
  if(MOD==0)pho<-par
  if(MOD==1)pho<-rep(par,m*(m-1)/2)
  pi<-Find.joint(p,pho,Omega = Omega,eps=eps)
  return(pi)
}

PredPi<-function(xt,piy,piz,h=1,Omega){
  Qmat<-matrix(NA,nrow(Omega),nrow(Omega))
  for(i in 1:nrow(Omega)){
  x0<-Omega[i,]
  if(h==1){
  i<-nrow(Omega)
  x0<-xt
  }
  for(j in 1:nrow(Omega)){
      x<-Omega[j,]
      indyx0<-which(x0==1)
      indzx0<-which(x0==0)
      if(length(indyx0)==0){ Qy<-1 }else{
        if(length(indyx0)==1){
          Qy<-sum(piy[Omega[,indyx0]==x[indyx0]])
        }else{
          match.idx <- which(apply(Omega[,indyx0], 1, all.equal, current=x[indyx0]) ==TRUE) 
          Qy<-sum(piy[match.idx]) 
        }}
        if(length(indzx0)==0){ Qz<-1 }else{
          if(length(indzx0)==1){
            Qz<-sum(piz[Omega[,indzx0]==x[indzx0]])
          }else{
            match.idx <- which(apply(Omega[,indzx0], 1, all.equal, current=x[indzx0]) ==TRUE) 
            Qz<-sum(piz[match.idx]) 
          }}
     Qmat[i,j]<-Qy*Qz 
    }
  }
  library(expm)  
  Qmath<-Qmat%^%h
  if(h==1){PredP<-Qmath[nrow(Omega),]}else{
    PredP<-Qmath[apply(Omega, 1, all.equal, current=xt)==TRUE,]}
  Pred<-Omega[which.max(PredP),]  
  return(list("Predicted.values"=Pred,"Predictive.Density"=PredP))
}

fbplot2<-function(fit,prob,mainplot="",ylim,xlab)
{
  fb<-fbplot(fit,plot=FALSE)
  ord.curves<-fit[,order(fb$depth,decreasing = TRUE)]
  med.curve<-fit[,which.max(fb$depth)]
  select.curves<-ord.curves[,1:ceiling(ncol(fit)*prob)]
  low<-apply(select.curves, 1, min)
  max<-apply(select.curves, 1, max)
  x<-0:(length(med.curve)-1)
  plot(x,med.curve,type="l",ylab = " ",xlab=xlab,lwd=2,
       xlim = c(0,8),main=mainplot,ylim = ylim)
  polygon(c(x,rev(x)), c(low,rev(max)), col = "deepskyblue", border = "gray", 
          lwd = 2)
  lines(x,med.curve,lwd=2)
}


