# Author: Ryan Yee
# Date: September 17, 2024
# Purpose: implementation of MALA to generate samples from a standard multivariate normal distribution
# Details: 
# Dependencies: 

### hyperparameters ###
lambda = 1
dim = 10
n_burn = 1000
n_post = 10000


### sampler ###

# randomly initialize chain
# x_init = rnorm(n = dim, mean = 0, sd = 1)
x_init = runif(n = dim, min = -5, max = 5)

x_samples = matrix(nrow = n_post, ncol = dim)
n_accept = 0

# MCMC loop
x = x_init
for (it in 1:(n_burn + n_post)){
  # proposal
  epsilon = rnorm(n = dim, mean = 0, sd = lambda)
  theta_star = (1 - (lambda^2 / 2) ) * x
  x_star = theta_star + epsilon
  
  # acceptance probability
  theta = (1 - (lambda^2 / 2) ) * x_star
  alpha_star = (prod(dnorm(x_star)) * prod(dnorm(x_star, mean = theta_star, sd = lambda))) / (prod(dnorm(x)) * prod(dnorm(x, mean = theta, sd = lambda)))
  alpha = min(1, alpha_star, na.rm = TRUE)
  
  # accept / reject
  if (alpha >= runif(1)){
    x = x_star
    if (it > n_burn) n_accept = n_accept + 1
  }
  
  # save sample
  if (it > n_burn) x_samples[it - n_burn, ] = x
}


### diagnostics ###
n_accept / n_post # acceptance rate
coda::effectiveSize(x_samples) # effective sample size
mean(x_samples)
sd(x_samples)

