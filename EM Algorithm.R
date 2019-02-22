## Fit the Location Mixture Model ##

library(MASS)
y<-galaxies/1000
N=length(y)
modes=c(4,6,8,11,15,20)
AIC=c()
BIC=c()

par(mfrow=c(2,3))
for(K in modes){
  
## Set initial values ##
  pim=c()
  mum=c()
  sm=var(y)
  pinikm<-matrix(nrow=K,ncol=N)
  for(k in 1:K){
    pim[k]=1/K
    mum[k]=min(y)+(k-1)*(max(y)-min(y))/(K-1)
  }
  
## Doing M Steps ##  
  for(i in 1:2000){
    
## Compute Pi_n,i,k ##
    for(n in 1:N){
      for(k in 1:K){
        pinikm[k,n]=pim[k]*dnorm(y[n],mum[k],sm^0.5)/(sum(pim*dnorm(y[n],mum,sm^0.5)))
      }
    }
## Update Pi and Mu ##
    for(k in 1:K){
      pim[k]=sum(pinikm[k,])/sum(pinikm)
      mum[k]=sum(pinikm[k,]*y)/sum(pinikm[k,])
    } 
## Update Sigma ##    
    sm=0
    for(n in 1:N){
      for(k in 1:K){
        sm=sm+pinikm[k,n]*(y[n]-mum[k])^2/sum(pinikm)
      }
    }
  }

## Graph Results ##  
  pdf<-function(x){
    y=0
    for(k in 1:K){
      y=y+pim[k]*dnorm(x,mum[k],sm^0.5)
    }
    return(y)
  }
  hist(y,breaks=30,freq=FALSE,main=paste("K=",K))
  curve(expr=pdf,col="blue",add=TRUE)
  
## Compute AIC & BIC ##
  l=0
  for(n in 1:N){
    l=l+log(sum(pim*dnorm(y[n],mum,sm^0.5)))
  }
  aic=2*(2*K+1)-2*l
  bic=(2*K+1)*log(N)-2*l
  AIC=c(AIC,aic)
  BIC=c(BIC,bic)
  
  print(paste("K =",K))
  print(pim)
  print(mum)
  print(sm)

}

AIC
BIC



## Fit the Location-Scale Mixture Model ##

library(MASS)
y<-galaxies/1000
N=length(y)
modes=c(3,4,5,6,7,8)
AIC=c()
BIC=c()

par(mfrow=c(2,3))
for(K in modes){
  
  ## Set initial values ##
  pim=c()
  mum=c()
  sm=c()
  pinikm<-matrix(nrow=K,ncol=N)
  for(k in 1:K){
    pim[k]=1/K
    mum[k]=min(y)+(k-1)*(max(y)-min(y))/(K-1)
    sm[k]=var(y)
  }
  
  ## Doing M Steps ##  
  for(i in 1:2000){
    
    ## Compute Pi_n,i,k ##
    for(n in 1:N){
      for(k in 1:K){
        pinikm[k,n]=pim[k]*dnorm(y[n],mum[k],sm[k]^0.5)/(sum(pim*dnorm(y[n],mum,sm^0.5)))
      }
    }
    ## Update Pi and Mu ##
    for(k in 1:K){
      pim[k]=sum(pinikm[k,])/sum(pinikm)
      mum[k]=sum(pinikm[k,]*y)/sum(pinikm[k,])
    }
    ## Update Sigma ##    
    for(k in 1:K){
      sm[k]=sum(pinikm[k,]*(y-mum[k])^2)/(sum(pinikm[k,]))
    }
  }
  
  ## Graph Results ##  
  pdf<-function(x){
    y=0
    for(k in 1:K){
      y=y+pim[k]*dnorm(x,mum[k],sm[k]^0.5)
    }
    return(y)
  }
  hist(y,breaks=30,freq=FALSE,main=paste("K=",K))
  curve(expr=pdf,col="blue",add=TRUE)
  
  ## Compute AIC & BIC ##
  l=0
  for(n in 1:N){
    l=l+log(sum(pim*dnorm(y[n],mum,sm^0.5)))
  }
  aic=2*(3*K)-2*l
  bic=3*K*log(N)-2*l
  AIC=c(AIC,aic)
  BIC=c(BIC,bic)
  
  print(paste("K =",K))
  print(pim)
  print(mum)
  print(sm)
}

AIC
BIC



## Fit location-scale mixtures of bivariate normals ##

library(datasets)
library(mvtnorm)
y<-faithful
N=272
modes=c(2,3,4,5)
AIC=c()
BIC=c()

par(mfrow=c(2,2))
for(K in modes){
  
  ## Set initial values ##
  pim<-c()
  mum<-matrix(nrow=K,ncol=2)
  sm<-array(dim=c(2,2,K))
  pinikm<-matrix(nrow=K,ncol=N)
  for(k in 1:K){
    pim[k]=1/K
    mum[k,]=c(min(y[,1])+(k-1)*(max(y[,1])-min(y[,1]))/(K-1), min(y[,2])+(k-1)*(max(y[,2])-min(y[,2]))/(K-1))
    sm[,,k]=var(y)
  }
  
  ## Doing M Steps ##  
  for(i in 1:200){
    
    ## Compute Pi_n,i,k ##
    for(n in 1:N){
      a=0
      for(j in 1:K){
        a=a+pim[j]*dmvnorm(y[n,],mum[j,],sm[,,j])
      }
      for(k in 1:K){
        pinikm[k,n]=pim[k]*dmvnorm(y[n,],mum[k,],sm[,,k])/a
      }
    }
    ## Update Pi ##
    for(k in 1:K){
      pim[k]=sum(pinikm[k,])/sum(pinikm)
    }
    ## Update Mu ##
    for(k in K){
      mum[k,1]=sum(pinikm[k,]*y[,1])/sum(pinikm[k,])
      mum[k,2]=sum(pinikm[k,]*y[,2])/sum(pinikm[k,])
    }
    ## Update Sigma ##    
    for(k in 1:K){
      sm[1,1,k]=sum(pinikm[k,]*(y[,1]-mum[k,1])^2)/sum(pinikm[k,])
      sm[2,2,k]=sum(pinikm[k,]*(y[,2]-mum[k,2])^2)/sum(pinikm[k,])
      sm[1,2,k]=sum(pinikm[k,]*(y[,1]-mum[k,1])*(y[,2]-mum[k,2]))/sum(pinikm[k,])
      sm[2,1,k]=sum(pinikm[k,]*(y[,1]-mum[k,1])*(y[,2]-mum[k,2]))/sum(pinikm[k,])
    }
  }
  
  ## Graph Results ##  
  plot(y,cex=0.5,main=paste("K=",K))
  pdf<-function(x1,x2){
    y=0
    for(k in 1:K){
      y=y+pim[k]*dmvnorm(c(x1,x2),mum[k,],sm[,,k])
    }
    return(y)
  }
  u=seq(40,100,0.1)
  v=seq(1.5,5.5,0.01)
  z<-outer(v,u,Vectorize(pdf))
  contour(v,u,z,add=TRUE)
  
  ## Compute AIC & BIC ##
  l=0
  for(n in 1:N){
    c=0
    for(k in 1:K){
      c=c+pim[k]*dmvnorm(y[n,],mum[k,],sm[,,k])
    }
    l=l+log(c)
  }
  aic=2*(6*K)-2*l
  bic=6*K*log(N)-2*l
  AIC=c(AIC,aic)
  BIC=c(BIC,bic)
  
  print(paste("K =",K))
  print(pim)
  print(mum)
  print(sm)
}

AIC
BIC


