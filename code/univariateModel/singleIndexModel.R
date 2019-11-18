# Function for analyzing results of a single acoustic index using the 
# model described in Doser et al (2019). 

# Author: Jeffrey W. Doser
# singleIndexModel.R ------------------------------------------------------

univar.sound.model <- function(n.iter, n.sites, data, ind, X, 
                               alpha.start, beta.start, sigma.sq.start, 
                               tau.sq.start, sim = FALSE, index = 'NDSI') {
  
  require(coda)
  
  # n.iter = number of iterations to run the Gibbs sampler
  # n.sites = number of recording site locations
  # data = the acoustic index data for each recording
  # ind = individual number for each soundscape recording
  # X = design matrix
  # alpha.start = starting value for all fixed regression coefficients
  # beta.start = starting value for random individual effects
  # sigma.sq.start = starting value for process error
  # tau.sq.start = starting value for random effect variance
  # index = acoustic index that is being used. Currently takes values 
  #         'H', 'AEI', 'NDSI', 'ACI', and 'sim' for simulated data.
  
  # Subroutines -----------------------------------------------------------
  # Multivariate normal random number generator
  rmvn <- function(n, mu=0, V = matrix(1)){
    p <- length(mu)
    if(any(is.na(match(dim(V),p))))
      stop("Dimension problem!")
    D <- chol(V)
    t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
  }
  
  # Inverse Gamma random number generator
  rigamma <- function(n, a, b){
    1/rgamma(n = n, shape = a, rate = b)
  } 
  
  # Logit transformation
  logit <- function(theta, a, b) {
    log((theta-a)/(b-theta))
  }
  
  # Inverse logit transformation
  logit.inv <- function(z, a, b) {
    b-(b-a)/(1+exp(z))
  }
  
  
  # Params ----------------------------------------------------------------
  
  # Will need to change accordingly when using different acoustic indices.
  y <- switch(index, 
              'NDSI' = as.matrix(logit(data, a = -1, b = 1)), 
              'H' = as.matrix(logit(data, a = 0, b = 1)), 
              'AEI' = as.matrix(logit(data, a = 0, b = 1)), 
              'ACI' = as.matrix(log(data)), 
              'sim' = as.matrix(data)
       )

  
  
  # For model validation, can remove different important parameters and see 
  # the effects this has on inference.
  
  # Remove the inherent differences parameter
  # X <- X[, -2]
  
  # Remove the rain parameter
  # X <- X[, -13]
  
  n.ind <- length(unique(ind))
  n <- nrow(X)
  W <- as.matrix(rep(1, n))
  n.alpha <- dim(X)[2]
  n.beta <- dim(W)[2]
  s <- rep(0, n.ind)
  for (i in 1:n.ind) {
    curr <- sum(ind == i)
    s[i] <- curr
  }
  
  # Priors ----------------------------------------------------------------
  
  # Normal non-informative prior on regression coefficients
  alpha.mu <- rep(0, n.alpha)
  alpha.var <- rep(10000, n.alpha)
  
  # Normal prior on random individual effects
  beta.mu <- matrix(0, nrow = n.ind, ncol = n.beta)
  
  # Inverse gamma prior on process variance
  sigma.sq.a <- 2
  sigma.sq.b <- 1
  
  # Inverse gamma prior on random effects variance
  tau.sq.a <- 2
  tau.sq.b <- 1
  
  # Starting values -------------------------------------------------------
  
  alpha <- rep(alpha.start, n.alpha)
  beta <- matrix(beta.start, nrow = n.ind, ncol = n.beta)
  sigma.sq <- sigma.sq.start
  tau.sq <- tau.sq.start
  
  
  # Sampler Prep ----------------------------------------------------------
  
  alpha.samples <- matrix(0, nrow = n.alpha, ncol = n.iter)
  beta.samples <- matrix(0, nrow = n.ind, ncol = n.iter)
  fitted.samples <- matrix(0, nrow = n, ncol = n.iter)
  tau.sq.samples <- rep(0, n.iter)
  sigma.sq.samples <- rep(0, n.iter)
  post.dens.samples <- matrix(0, n, n.iter)
  post.log.dens.samples <- matrix(0, n, n.iter)
  XTX <- t(X)%*%X
  
  # Gibbs Sampler ---------------------------------------------------------
  
  for (k in 1:n.iter) {
    
    # Update alpha samples
    V <- chol2inv(chol(XTX / sigma.sq + diag(1/alpha.var, n.alpha)))
    sum.diffs <- rep(0, n.alpha)
    for (i in 1:n.ind) {
      X.i <- X[ind == i, ]
      y.i <- y[ind == i]
      # W.i <- W[ind == i, ]
      W.i <- as.matrix(W[ind == i, ])
      # sum.diffs <- sum.diffs + t(X.i) %*% (y.i - ones.s*beta[i] - w.i * beta.1[i])
      sum.diffs <- sum.diffs + t(X.i) %*% (y.i - W.i%*% beta[i, ])
    }
    v <- sum.diffs / sigma.sq
    alpha <- rmvn(1, V%*%v, V)
    
    # Update beta samples
    for ( i in 1:n.ind) {
      X.i <- X[ind == i, ]
      y.i <- y[ind == i]
      # W.i <- W[ind == i, ]
      W.i <- as.matrix(W[ind == i, ])
      V <- chol2inv(chol(t(W.i)%*%W.i / sigma.sq + diag(1/tau.sq, n.beta)))
      v <- (t(W.i) %*% (y.i - X.i %*% alpha)) / sigma.sq
      beta[i, ] <- rmvn(1, V%*%v, V)
    }
    
    # Update tau.sq samples
    a <- tau.sq.a + .5 * n.ind * n.beta
    b <- tau.sq.b + .5 * sum(t(beta)%*%beta)
    tau.sq <- rigamma(1, a, b)
    
    # Update sigma.sq samples
    a <- sigma.sq.a + .5 * sum(s)
    b.sum <- 0
    for (i in 1:n.ind) {
      X.i <- X[ind == i, ]
      y.i <- y[ind == i]
      # W.i <- W[ind == i, ]
      W.i <- as.matrix(W[ind == i, ])
      b.sum <- b.sum + t((y.i - (X.i%*%alpha + W.i%*%beta[i, ]))) %*% 
        (y.i - (X.i%*%alpha + W.i%*%beta[i, ]))
    }
    b <- sigma.sq.b + .5 * b.sum
    sigma.sq <- rigamma(1, a, b)
    
    # Get values for WAIC -------------------------------------------------
    mu <- X %*% alpha + beta[ind]
    post.dens.samples[, k] <- dnorm(y, mu, sigma.sq)
    post.log.dens.samples[, k] <- dnorm(y, mu, sigma.sq, log = TRUE)
    
    # Save Samples
    alpha.samples[, k] <- alpha
    beta.samples[, k] <- beta[, 1]
    tau.sq.samples[k] <- tau.sq
    sigma.sq.samples[k] <- sigma.sq
    
    print(paste(k/n.iter * 100, " percent complete", sep = ""))
    
  }
  
  # Output ----------------------------------------------------------------
  return(
    list(
      alpha.samples = alpha.samples, beta.samples = beta.samples, 
      tau.sq.samples = tau.sq.samples, sigma.sq.samples = sigma.sq.samples,
      post.dens.samples = post.dens.samples, 
      post.log.dens.samples = post.log.dens.samples, n.iter = n.iter
    )
  )
}
