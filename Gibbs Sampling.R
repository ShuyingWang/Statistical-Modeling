## Importance Sampling ##
# (a) #
par(mfrow=c(1,2))
n=100
z=rnorm(10000,0,3^0.5)  # sampling from normal distribution
w=dt(z,3)/dnorm(z,0,3^0.5)/sum(dt(z,3)/dnorm(z,0,3^0.5))  # compute the weights for resampling
theta=sample(z,n,replace=FALSE,prob=w)  # resampling from z_i's
hist(theta,breaks=20,freq=FALSE,main="n=100")
curve(dt(x,3),add=TRUE)

n=10000
z=rnorm(100000,0,3^0.5)  # sampling from normal
w=dt(z,3)/dnorm(z,0,3^0.5)/sum(dt(z,3)/dnorm(z,0,3^0.5))  # compute the weights for resampling
theta=sample(z,n,replace=FALSE,prob=w)  # resampling from z_i's
hist(theta,breaks=50,freq=FALSE,main="n=10,000")
curve(dt(x,3),add=TRUE)

# (b) #
n=100 # repeat for n=10000
theta=rnorm(n,0,3^0.5) # sampling from normal distribution
g=theta*dt(theta,3)/dnorm(theta,0,3^0.5)
h=theta^2*dt(theta,3)/dnorm(theta,0,3^0.5)
Etheta=mean(g)  # estimate expectation of theta
Vtheta=mean(h)-Etheta^2  # estimate variance of theta
Etheta
Vtheta


## Gibbs Sampling ##
library(mvtnorm)
th=c(0,2) # set the initial value theta_0
sigma=0.4375^0.5

# The first 1000 iterations before burn in
for(i in 1:1000){
  th[1]=rnorm(1,0.75*th[2]-1.5,sigma) # sample theta_1 from p(theta_1 | theta_2)
  th[2]=rnorm(1,0.75*th[1]+2,sigma) # sample theta_2 from p(theta_2 | theta_1)
}
theta=matrix(nrow=1000,ncol=2)

# 1000 thinned sample after burn in
for(i in 1:1000){
  for(j in 1:5){
    th[1]=rnorm(1,0.75*th[2]-1.5,sigma) # sample theta_1 from p(theta_1 | theta_2)
    th[2]=rnorm(1,0.75*th[1]+2,sigma) # sample theta_2 from p(theta_2 | theta_1)
  }
  theta[i,]=th
}

# Graph scatter plot of the sample, superimposed over a contour of the actual density
plot(theta,cex=0.5)
s=matrix(c(1,0.75,0.75,1),nrow=2,ncol=2)
pdf<-function(x1,x2){
  y=dmvnorm(c(x1,x2),c(0,2),s)
  return(y)
}
u=seq(-1.5,5.5,0.1)
v=seq(-3.5,3.5,0.1)
z<-outer(v,u,Vectorize(pdf))
contour(v,u,z,add=TRUE)

# Graph the trace and autocorrelation plot of sample #
par(mfrow=c(2,2))
plot(theta[,1],type="l",main="theta_1")
acf(theta[,1],main="theta_1")
plot(theta[,2],type="l",main="theta_2")
acf(theta[,2],main="theta_2")


## AR(1) Process ##
y=c()
y[1]=rnorm(1,0,1)  # sample the first data from N(0,1)

# Simulate AR(1) process
for(t in 1:499){
  y[t+1]=0.75*y[t]+rnorm(1,0,1) 
}
# Compute the posterior
a=sum(y^2)-y[1]^2
b=sum(y^2)-y[500]^2
c=0
for(t in 1:499){
  c=c+y[t]*y[t+1]
}
mu=c/(b+1)
nu=500
k=b+1
s=(2+a-c^2/(b+1))/500

# Gibbs sampler from posterior
library(invgamma)
rho=c()
sigma=c()
r=mu
sig=1 # set the initial value
# The first 1000 iterations before burn in
for(i in 1:1000){
  r=rnorm(1,mu,sig/k) 
  sig=rinvgamma(1,(nu+1)/2,(k*(r-mu)^2+nu*s)/2) 
}
# 2000 thinned sample after burn in
for(i in 1:2000){
  for(j in 1:5){
    r=rnorm(1,mu,sig/k) 
    sig=rinvgamma(1,(nu+1)/2,(k*(r-mu)^2+nu*s)/2) 
  }
  rho=c(rho,r)
  sigma=c(sigma,sig)
}
# Credible interval
quantile(rho,probs=0.05)
quantile(rho,probs=0.95)
quantile(sigma,probs=0.05)
quantile(sigma,probs=0.95)


## Location-scale mixture normal ##
library(MASS)
library(extraDistr)
library(DirichletReg)
library(invgamma)
library(graphics)
y<-galaxies/1000
N=length(y)
modes=c(3,4,5,6,7,8)
# Choose prior hyper-parameters 
mu0=mean(y)
s0=var(y)/2
a0=1
b0=1
par(mfrow=c(2,3))
for(K in modes){
  # Gibbs sampler
  mu=c()
  s=c()
  pi=c()
  n=c()
  alpha=c()
  piz=c()
  B=100
  z=matrix(nrow=N,ncol=K)
  pi_s=matrix(nrow=B,ncol=K)
  mu_s=matrix(nrow=B,ncol=K)
  s_s=matrix(nrow=B,ncol=K)
  
  # Set initial values
  for(k in 1:K){
    mu[k]=min(y)+(k-1)*(max(y)-min(y))/(K-1)
    s[k]=var(y)
    pi[k]=1/K
  } 
  
  # First 1000 iterations before burn in
  for(t in 1:100){
    # sample z_i's
    for(i in 1:N){
      piz=(pi*dnorm(y[i],mu,s^0.5))/sum(pi*dnorm(y[i],mu,s^0.5))
      z[i,]=rmnom(1, 1, piz)
    } 
      for(k in 1:K){
        n[k]=sum(z[,k])
        alpha[k]=n[k]+1
      }
    # sample pi
    pi=rdirichlet(1,alpha) 
    # sample mu_k's and sigma_k's
    for(k in 1:K){
      skn=1/(1/s0+n[k]/s[k])
      mkn=skn*(mu0/s0+sum(y*z[,k])/s[k])
      mu[k]=rnorm(1,mkn,skn^0.5)
      a=a0+n[k]/2
      b=b0+(sum((y-mu[k])^2*z[,k]))/2
      s[k]=rinvgamma(1,a,b)
    }
  }
  # samples after burn in
  for(t in 1:B){
      for(i in 1:N){
        piz=(pi*dnorm(y[i],mu,s^0.5))/sum(pi*dnorm(y[i],mu,s^0.5))
        z[i,]=rmnom(1, 1, piz)
      } 
      for(k in 1:K){
        n[k]=sum(z[,k])
        alpha[k]=n[k]+1
      }
      pi=rdirichlet(1,alpha) 
      for(k in 1:K){
        skn=1/(1/s0+n[k]/s[k])
        mkn=skn*(mu0/s0+sum(y*z[,k])/s[k])
        mu[k]=rnorm(1,mkn,skn^0.5)
        a=a0+n[k]/2
        b=b0+(sum((y-mu[k])^2*z[,k]))/2
        s[k]=rinvgamma(1,a,b)
      }
    mu_s[t,]=mu
    pi_s[t,]=pi
    s_s[t,]=s
  }
  
  # compute mean density and credible interval
  hist(y,breaks=30,freq=FALSE,main=paste("K=",K))
  u=seq(5,40,0.01)
  nu=length(u)
  w=matrix(nrow=nu,ncol=B)
  lb=c()
  ub=c()
  mead=c()
  for(i in 1:nu){
    for(j in 1:B){
      v=0
      for(k in 1:K){
        v=v+pi_s[j,k]*dnorm(u[i],mu_s[j,k],s_s[j,k]^0.5)
      }
      w[i,j]=v
    }
    lb[i]=quantile(w[i,],probs=0.05)
    ub[i]=quantile(w[i,],probs=0.95)
    mead[i]=mean(w[i,])
  }
  lines(x=u,y=mead,col="blue")
  polygon(x=c(u,rev(u)),y=c(lb,rev(ub)),col=rgb(0.1,0.1,0.1,0.2),border=NA)
}  



## Location-scale mixture bivariate normal ##

library(datasets)
library(mvtnorm)
library(extraDistr)
library(DirichletReg)
library(LaplacesDemon)
y<-faithful
N=272
y=array(c(y[,1],y[,2]),c(N,2))
modes=c(2,3,4,5)
# choose prior hyper-parameters
mu0=c(mean(y[,1]),mean(y[,2]))
s0=var(y)
nu0=2
psi0=matrix(c(2,0,0,2),nrow=2,ncol=2)
B=200

par(mfrow=c(2,2))
for(K in modes){
  # Set initial values
  pi=c()
  mu=matrix(nrow=K,ncol=2)
  s<-array(dim=c(2,2,K))
  n=c()
  alpha=c()
  piz=c()
  z=matrix(nrow=N,ncol=K)
  S=matrix(c(0,0,0,0),nrow=2,ncol=2)
  
  for(k in 1:K){
    pi[k]=1/K
    mu[k,]=c(min(y[,1])+(k-1)*(max(y[,1])-min(y[,1]))/(K-1), min(y[,2])+(k-1)*(max(y[,2])-min(y[,2]))/(K-1))
    s[,,k]=var(y)
  }
  # 1000 iterations before burn
  for(t in 1:1000){
    # sample z_i's
    for(i in 1:N){
      for(k in 1:K){
        piz[k]=pi[k]*dmvnorm(y[i,],mu[k,],s[,,k])
      }
      piz=piz/sum(piz)
      z[i,]=rmnom(1, 1, piz)
    } 
    for(k in 1:K){
      n[k]=sum(z[,k])
      alpha[k]=n[k]+1
    }
    # sample pi
    pi=rdirichlet(1,alpha) 
    # sample mu_k's and sigma_k's
    for(k in 1:K){
      skn=solve(solve(s0)+n[k]*solve(s[,,k]))
      mkn=skn%*%(solve(s0)%*%mu0+solve(s[,,k])%*%t(y)%*%z[,k])
      mu[k,]=rmvnorm(1,mkn,skn)
      S=matrix(c(0,0,0,0),nrow=2,ncol=2)
      for(i in 1:N){
        S=S+z[i,k]*(y[i,]-mu[k,])%*%t(y[i,]-mu[k,])
      }
      s[,,k]=rinvwishart(nu0+n[k], S+psi0)
    }
  }
  # Density function
  pdf<-function(x1,x2){
    y=0
    for(k in 1:K){
      y=y+pi[k]*dmvnorm(c(x1,x2),mu[k,],s[,,k])
    }
    return(y)
  }
  u=seq(1.5,5.5,0.1)
  v=seq(40,100,1)
  un=length(u)
  vn=length(v)
  w=array(dim=c(un,vn,B))
  
  # sample after burn in
  for(t in 1:B){
    # sample z_i's
    for(i in 1:N){
      for(k in 1:K){
        piz[k]=pi[k]*dmvnorm(y[i,],mu[k,],s[,,k])
      }
      piz=piz/sum(piz)
      z[i,]=rmnom(1, 1, piz)
    } 
    for(k in 1:K){
      n[k]=sum(z[,k])
      alpha[k]=n[k]+1
    }
    # sample pi
    pi=rdirichlet(1,alpha) 
    # sample mu_k's and sigma_k's
    for(k in 1:K){
      skn=solve(solve(s0)+n[k]*solve(s[,,k]))
      mkn=skn%*%(solve(s0)%*%mu0+solve(s[,,k])%*%t(y)%*%z[,k])
      mu[k,]=rmvnorm(1,mkn,skn)
      S=matrix(c(0,0,0,0),nrow=2,ncol=2)
      for(i in 1:N){
        S=S+z[i,k]*(y[i,]-mu[k,])%*%t(y[i,]-mu[k,])
      }
      s[,,k]=rinvwishart(nu0+n[k], S+psi0)
    }
    for(i in 1:un){
      for(j in 1:vn){
        w[i,j,t]=pdf(u[i],v[j])
      }
    }
  }
  # Compute and graph mean density
  h=matrix(nrow=un,ncol=vn)
  for(i in 1:un){
    for(j in 1:vn){
      h[i,j]=mean(w[i,j,])
    }
  }
  plot(y,cex=0.5,main=paste("K=",K))
  contour(u,v,h,add=TRUE)
}