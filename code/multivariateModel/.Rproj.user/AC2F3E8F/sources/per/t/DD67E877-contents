# This file contains the function multivar.sound.model that runs the 
# multivariate Gibbs sampler described in Doser et al (2019). 

# Author: Jeffrey W. Doser

multivar.sound.model <- function(n.iter, fileName, geoFileName, indexNames, 
                                 yearsPostTrt, alpha.start, beta.start, 
                                 sigma.sq.start, lambda.start, sim = FALSE) {
  
  require(coda)
  require(dplyr)
  source("computeDesignMatrix.R")
  # Multivariate normal random number generator
  rmvn <- function(n, mu=0, V = matrix(1)) {
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
  
  # Logit data transformation
  logit <- function(theta, a, b) {
    log((theta-a)/(b-theta))
  }
  
  # Inverse logit transformation
  logit.inv <- function(z, a, b) {
    b-(b-a)/(1+exp(z))
  }
  
  
  # Params ------------------------------------------------------------------
  out <- computeDesignMatrix(fileName, geoFileName, 
                             indexNames = indexNames, 
                             yearsPostTrt = yearsPostTrt)
  X <- out$X
  dat <- out$dat
  n.sites <- length(unique(dat$recording.site))
  indices <- dat$value
  y <- as.matrix(dat$value)
  n.indices <- out$n.indices
  ind <- dat$ind
  n.ind <- length(unique(ind))
  n <- nrow(X)
  n.alpha <- dim(X)[2]
  n.beta <- 1
  s <- rep(0, n.ind)
  for (i in 1:n.ind) {
    curr <- sum(ind == i)
    s[i] <- curr / n.indices
  }
  ind.beta <- rep(1:n.ind, each = n.indices)
  
  # Priors ------------------------------------------------------------------
  
  # Normal prior on fixed regression coefficients
  alpha.mu <- rep(0, n.alpha)
  alpha.var <- rep(10000, n.alpha)
  
  # Normal prior on random individual effects
  beta.mu <- matrix(0, nrow = n.ind * n.indices, ncol = n.beta)
  
  # Inverse gamma prior on process variance
  sigma.sq.a <- 2
  sigma.sq.b <- 1
  
  # Inverse wishart prior on 10 x 10 covariance matrix
  R <- diag(0.1, n.indices)
  r <- n.indices
  
  
  # Starting values ---------------------------------------------------------
  alpha <- rep(alpha.start, n.alpha)
  beta <- as.matrix(rep(beta.start, n.ind * n.indices))
  sigma.sq <- sigma.sq.start
  lambda <- diag(lambda.start, n.indices)
  y.hat <- rep(0, n)
  
  
  # Sampler Prep ------------------------------------------------------------
  
  alpha.samples <- matrix(0, nrow = n.alpha, ncol = n.iter)
  beta.samples <- matrix(0, nrow = n.ind * n.indices, ncol = n.iter)
  fitted.samples <- matrix(0, nrow = n, ncol = n.iter)
  lambda.samples <- matrix(0, nrow = n.indices * n.indices, ncol = n.iter)
  sigma.sq.samples <- rep(0, n.iter)
  XTX <- t(X)%*%X
  
  # Gibbs Sampler -----------------------------------------------------------
  
  for (k in 1:n.iter) {
    
    # Update alpha samples
    V <- chol2inv(chol(XTX / sigma.sq + diag(1/alpha.var, n.alpha)))
    sum.diffs <- rep(0, n.alpha)
    for (i in 1:n.ind) {
      X.i <- X[ind == i, ]
      y.i <- y[ind == i]
      # Faster way to code the indicator matrix described in the sampler
      beta.i <- rep(beta[ind.beta ==i], times = s[i]) 
      sum.diffs <- sum.diffs + t(X.i) %*% (y.i - beta.i)
    }
    v <- sum.diffs / sigma.sq
    alpha <- rmvn(1, V%*%v, V)
    
    # Update beta samples
    for (i in 1:n.ind) {
      X.i <- X[ind == i, ]
      y.i <- y[ind == i]
      V <- chol2inv(chol(diag(s[i] / sigma.sq, n.indices) + chol2inv(chol(lambda))))
      # Indicator matrix (aka w)
      ones.mat <- diag(1, n.indices)
      ones.mat <- matrix(rep(ones.mat, s[i]), nrow = n.indices * s[i], 
                         ncol = n.indices, byrow = T)
      v <-  t(ones.mat) %*% (y.i - X.i %*% alpha) / sigma.sq
      beta[((i-1)*n.indices + 1):(i*n.indices), ] <- rmvn(1, V%*%v, V)
    }
    
    # Update sigma.sq samples
    a <- sigma.sq.a + 5 * sum(s)
    b.sum <- 0
    for (i in 1:n.ind) {
      X.i <- X[ind == i, ]
      y.i <- y[ind == i]
      # This is a faster way to code the identicator matrix described in the sampler
      beta.i <- rep(beta[ind.beta ==i], times = s[i]) 
      b.sum <- b.sum + t((y.i - (X.i%*%alpha + beta.i))) %*% 
        (y.i - (X.i%*%alpha + beta.i))
    }
    b <- sigma.sq.b + .5 * b.sum
    sigma.sq <- rigamma(1, a, b)
    
    # Update lambda (inverse) samples 
    a <- n.ind + r
    beta.sum <- matrix(0, ncol = n.indices, nrow = n.indices)
    for (i in 1:n.ind) {
      beta.sum <- beta.sum + (beta[((i-1)*n.indices + 1):(i*n.indices), ] %*% 
                                t(beta[((i-1)*n.indices + 1):(i*n.indices), ]))
    }
    b <- chol2inv(chol(beta.sum + r*R))
    # Note that the rWishart function returns a 3D array, so you need to 
    # extract only the first two dimensions
    lambda <- chol2inv(chol(rWishart(1, a, b)[, , 1]))
    
    # Save Samples
    alpha.samples[, k] <- alpha
    beta.samples[, k] <- beta[, 1]
    lambda.samples[, k] <- c(lambda)
    sigma.sq.samples[k] <- sigma.sq
    
    print(paste(k/n.iter * 100, " percent complete", sep = ""))
    
  }
  
  # Output
  return(
    list(
      alpha.samples = alpha.samples, beta.samples = beta.samples, 
      lambda.samples = lambda.samples, sigma.sq.samples = sigma.sq.samples
    )
  )
  
  
}

