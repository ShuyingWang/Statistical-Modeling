
### Forward-Backward Algorithm and Gibbs sampling for Hidden Markov Models 

### Bayesian HMM with location-scale Normal mixtures as emission densities

library(datasets)
library(MASS)
library(extraDistr)
library(DirichletReg)
library(invgamma)
library(graphics)
library(RDS)

y <- faithful[ ,2]  
N <- length(y)
B <- 500  # Number of simulations
modes <- c(2, 3)  # Number of components in Normal mixture model



## Set prior hyper-parameters -----------------------------------------------------------------------------------------------------

mu_0 <- mean(y)
sigma_0 <- var(y)/2  # Hyper-parameters for mu
a_0 <- 1
b_0 <- 1  # Hyper-parameters for sigma


par(mfrow = c(2, 2))

for (K in modes){
  
## Set initial values for parameters ---------------------------------------------------------------------------------------------------------
  
  pi_0 <- replicate(K, 1/K)  # Initial distribution for latent variable z -- p(z_1)
  P <- matrix(1/K, nrow = K, ncol = K)  # Transition ptobability matrix for z -- p(z_t | z_(t-1))
  mu <- c()
  for (k in 1:K){  
    mu[k] = min(y) + (k-1) * (max(y) - min(y)) / (K-1)
  }  
  sigma <- replicate(K, var(y))  # Emission distibution for y given z -- p(y_t | z_t)


    
## Define function for forward-backward algorith (FBA) ---------------------------------------------------------------------------
  
  FBA <- function(){
    
    # compute normalized forward messages
    
    alpha <- matrix(nrow=N, ncol=K)
    alpha[1, ] = dnorm(y[1], mu, sigma^0.5) * pi_0
    alpha[1, ] = alpha[1, ] / sum(alpha[1, ]) 
    
    for (t in 2:N){
      for (k in 1:K){
        alpha[t, k] = dnorm(y[t], mu[k], sigma[k]^0.5) * sum(P[ ,k] * alpha[t-1, ])
      }
      alpha[t, ] = alpha[t, ] / sum(alpha[t, ])
    }
    
    # Sample backwards from z_T to z_1
    
    z <- c()  
    Z <- matrix(nrow=N, ncol=K)  # Latent variable z
    p_zt <- c()  # Conditional probability for sampling z backwards -- p(z_t | z_(t+1):T)
    
    p_zt = alpha[N, ] / sum(alpha[N, ])  
    Z[N,] = rmnom(1, 1, p_zt)
    z[N] = Z[N,] %*% c(1:K)  
    
    for (t in (N-1):1){
      p_zt = alpha[t, ] * P[ ,z[t+1]]
      p_zt = p_zt / sum(p_zt)
      Z[t, ] = rmnom(1, 1, p_zt)
      z[t] = Z[t, ] %*% c(1:K)
    }  
    return(Z)
  }  
  

    
## Use Gibbs sampling to sample parameters from posterior distribution ----------------------------------------------------------
  
  # Set up for Gibbs sampling
  
  pi_s <- matrix(nrow=B, ncol=K)  
  mu_s <- matrix(nrow=B, ncol=K)
  sigma_s <- matrix(nrow=B, ncol=K)  # Matrices to contain stationary samples of parameters
  pi_predict <- matrix(nrow=B, ncol=K) # Matrix to contain the predictive density of latent variable
  
  
  # The first 1000 iterations before burn in
  for (i in 1:1000){
    
    # sample z_1:T
  
    Z = FBA()
    z = Z %*% c(1:K)
    n_k <- apply(Z, 2, sum)  # Count number of z_t = k 
    
    n_jk <- matrix(0, nrow = K, ncol = K)
    for (t in 2:N){
      n_jk[z[t-1], z[t]] = n_jk[z[t-1], z[t]] + 1  # Count number of pairs that z_(t-1) = j, z_t = k
    }
    
    # sample mu_k's and sigma_k's
    
    for (k in 1:K){
      v = 1 / (1 / sigma_0 + n_k[k] / sigma[k]) 
      m = v * (mu_0 / sigma_0 + sum(y*Z[ ,k]) / sigma[k])  # mean and variance for the conditional distribution of mu_k
      mu[k] = rnorm(1, m, v^0.5)  
      
      a = a_0 + n_k[k]/2
      b = b_0 + (sum((y-mu[k])^2 * Z[,k])) * 0.5  # parameters for the conditional distribution of sigma_k
      sigma[k] = rinvgamma(1, a, b)
    }
    
    # sample pi_0 and pi_j
  
    pi_0 = rdirichlet(1, replicate(K, 1/K) + Z[1, ])
    
    for (j in 1:K){
      P[j, ] = rdirichlet(1, replicate(K, 1/K) + n_jk[j, ]) 
    }
  }
  
  
  # Continue to sample after burn in
  
  for (i in 1:B){

    # sample z_1:T
    
    Z = FBA()
    z = Z %*% c(1:K)
    n_k <- apply(Z, 2, sum)  # Count number of z_t = k 
    
    n_jk <- matrix(0, nrow = K, ncol = K)
    for (t in 2:N){
      n_jk[z[t-1], z[t]] = n_jk[z[t-1], z[t]] + 1  # Count number of pairs that z_(t-1) = j, z_t = k
    }
    
    # sample mu_k's and sigma_k's
    
    for (k in 1:K){
      v = 1 / (1 / sigma_0 + n_k[k] / sigma[k]) 
      m = v * (mu_0 / sigma_0 + sum(y*Z[ ,k]) / sigma[k])  # mean and variance for the conditional distribution of mu_k
      mu[k] = rnorm(1, m, v^0.5)  
      
      a = a_0 + n_k[k]/2
      b = b_0 + (sum((y-mu[k])^2 * Z[,k])) * 0.5  # parameters for the conditional distribution of sigma_k
      sigma[k] = rinvgamma(1, a, b)
    }
    
    # sample pi_0 and pi_j
    
    pi_0 = rdirichlet(1, replicate(K, 1/K) + Z[1, ])
    
    for (j in 1:K){
      P[j, ] = rdirichlet(1, replicate(K, 1/K) + n_jk[j, ]) 
    }

    # record samples in the current iteration
    mu_s[i, ] = mu
    pi_s[i, ] = get.stationary.distribution(P)
    sigma_s[i, ] = sigma
    pi_predict[i, ] = P[z[N], ]
  }
  
  
  
## Compute stationary & predictive density and credible interval ---------------------------------------------------------------------------------------
  
  u <- seq(40, 100, 0.01)
  n <- length(u)
  sample_densities <- matrix(nrow = n, ncol=B)
  low_bound <- c() 
  up_bound <- c()
  station_density <- c()
  predict_density <- c()
  
  # compute stationary density
  for (i in 1:n){
    for (j in 1:B){
      sample_densities[i, j] = sum(pi_s[j, ] * dnorm(u[i], mu_s[j, ], sigma_s[j, ]^0.5))
    }
    low_bound[i] = quantile(sample_densities[i, ], probs = 0.05)  
    up_bound[i] = quantile(sample_densities[i, ], probs = 0.95)  # Compute upper and lower bound of the credible interval
    station_density[i] = mean(sample_densities[i, ])  # Compute the mean stationary density
  }
  
  # Graph the histogram of data, superimposed with the stationary density and credible interval
  
  hist(y, breaks = 30, freq = FALSE, main = paste("K = ", K), ylab = ("stationary density"))
  lines(x = u, y = station_density, col = "blue")
  polygon(x = c(u, rev(u)), y = c(low_bound, rev(up_bound)), col = rgb(0.1, 0.1, 0.1, 0.2), border = NA)
  
  
  
  # Compute predictive densities 
  for (i in 1:n){
    for (j in 1:B){
      sample_densities[i, j] = sum(pi_predict[j, ] * dnorm(u[i], mu_s[j, ], sigma_s[j, ]^0.5))
    }
    low_bound[i] = quantile(sample_densities[i, ], probs = 0.05)
    up_bound[i] = quantile(sample_densities[i, ],  probs= 0.95)
    predict_density[i] = mean(sample_densities[i, ])
  }
  
  # Graph the histogram of data, superimposed with the predictive density and credible interval
  
  hist(y, breaks = 30, freq = FALSE, main = paste("K=", K), ylab = ("predictive density"))
  lines(x = u, y = predict_density, col="blue")
  polygon(x = c(u, rev(u)), y = c(low_bound, rev(up_bound)), col = rgb(0.1, 0.1, 0.1, 0.2), border=NA)
}  