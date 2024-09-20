# Author: Ryan Yee
# Date: September 19, 2024
# Purpose: functions to test different settings for MALA standard multivariate normal sampler
# Details: 
# Dependencies: mcmcse

# inputs:
#   dim - dimension of samples
#   lambda - MALA step size
#   n_burn - number of burn-in iterations
#   n_post - number of posterior samples returned
#   save_samples - TRUE if you want the function to return the actual samples
# output: list with the following attributes
#   matrix of size (n_post x dim) of posterior samples from MALA sampler is save_samples = TRUE
#   sample mean
#   sample variance
#   effective sample size
#   acceptance rate
mala_std_mvn_sample = function(dim, lambda, n_burn, n_post, save_samples = FALSE){
  
  # output containers
  samples = matrix(nrow = n_post, ncol = dim)
  n_accept = 0
  
  # randomly initialize chain
  x = rnorm(n = dim, mean = 0, sd = 1)
  
  # MCMC loop
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
    if (it > n_burn) samples[it - n_burn, ] = x
  }
  
  # diagnostics
  diag = list(
    mean_bias = sqrt(sum(apply(X = samples, MARGIN = 2, FUN = mean)^2)),
    var_bias = sqrt(sum((apply(X = samples, MARGIN = 2, FUN = sd)^2 - rep(1, times = dim))^2)),
    ess = mcmcse::multiESS(samples),
    acceptance_rate = n_accept / n_post
  )
  
  if (save_samples){
    out = list(diag, samples = samples)
  } else {
    out = diag
  }
  
  return(out)
}

