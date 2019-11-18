# Function for simulating data for testing the multivariate model 
# described in Doser et al (2019)

# Author: Jeffrey W. Doser

# Multivariate normal random number generator -----------------------------
sim.multi.dat <- function(n.sites = 13, n.dc = 4, n.days = 30, 
                          n.years = 5, n.indices = 10, n.trt = 4) {
  
  rmvn <- function(n, mu=0, V = matrix(1)){
    p <- length(mu)
    if(any(is.na(match(dim(V),p))))
      stop("Dimension problem!")
    D <- chol(V)
    t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
  }
  
  # Simulate Data -----------------------------------------------------------
  
  set.seed(1)
  n <- n.sites * n.dc * n.days * n.years * n.indices
  n.ind <- n.days * n.dc * n.sites * n.indices
  n.ind.one <- n.days * n.dc * n.sites
  
  n.alpha <- ((n.years - 1) * 2 + 2) * n.indices
  
  # Fixed effect coefficients
  alpha <- rnorm(n.alpha, 0, 2)
  
  # Process variance
  sigma.sq <- 1
  Sigma.sq <- diag(sigma.sq, n.indices)
  
  # Random effect covariance matrix
  # Can manually change these values as desired
  lambda <- matrix(0, nrow = n.indices, ncol = n.indices)
  diag(lambda) <- 1.5
  lambda[2, 2] <- .4
  lambda[3, 3] <- .9
  # Individual random effects
  n.beta <- 1
  beta.mu <- rep(0, n.indices)
  beta.mat <- rmvn(n.ind.one, beta.mu, lambda)
  beta <- c(beta.mat)
  # beta.1 <- rnorm(n.ind, 0, sd = sqrt(tau.1))
  # beta <- cbind(beta.0, beta.1)
  
  n.ctrl <- n.sites - n.trt
  
  
  # Construct design matrix -------------------------------------------------
  
  # Create treatment groups
  trt.ind <- rep(0, n)
  trt.ind[1:(6000*n.trt)] <- 1
  year <- rep(rep(c(2013, 2014, 2015, 2016, 2017), each = n.indices), length.out = n)
  index.ind <- rep(1:n.indices, length.out = n)
  X <- matrix(0, ncol = ((n.years - 1) * 2 + 2) * n.indices, nrow = n)
  num <- ncol(X) / n.indices
  for (i in 1:n.indices) {
    # Respective intercept for each 
    X[, (i-1)*num + 1] <- ifelse(index.ind == i, 1, 0)
    X[, (i-1)*num + 2] <- ifelse(index.ind == i & trt.ind == 1, 1, 0)
    X[, (i-1)*num + 3] <- ifelse(index.ind == i & year == 2014, 1, 0)
    X[, (i-1)*num + 4] <- ifelse(index.ind == i & year == 2015, 1, 0)
    X[, (i-1)*num + 5] <- ifelse(index.ind == i & year == 2016, 1, 0)
    X[, (i-1)*num + 6] <- ifelse(index.ind == i & year == 2017, 1, 0)
    X[, (i-1)*num + 7] <- ifelse(index.ind == i & year == 2014 & trt.ind == 1, 1, 0)
    X[, (i-1)*num + 8] <- ifelse(index.ind == i & year == 2015 & trt.ind == 1, 1, 0)
    X[, (i-1)*num + 9] <- ifelse(index.ind == i & year == 2016 & trt.ind == 1, 1, 0)
    X[, (i-1)*num + 10] <- ifelse(index.ind == i & year == 2017 & trt.ind == 1, 1, 0)
  }
  
  # Indicator for individual
  ind <- rep(1:n.ind.one, each = n.years * n.indices)
  
  # Construct likelihood
  I.s <- diag(1, n.years)
  V <- kronecker(I.s, Sigma.sq)
  y <- rep(0, n)
  for (i in 1:n.ind.one) {
    x.i <- X[ind == i, ]
    # W.i <- W[ind == i, ]
    # v <- x.i %*% alpha + rep(1, n.years)*beta.0[i] + w.i * beta.1[i]
    v <- x.i %*% alpha + beta[((i-1)*n.indices + 1):(((i-1)*n.indices) + n.indices)]
    y[((i-1)*n.years*n.indices + 1):(((i-1)*n.years*n.indices)+n.years*n.indices)] <- rmvn(n = 1, mu = v, V = V)
  }
  
  return(
    list(
      X = X, y = y, beta = beta, alpha = alpha, n.sites = n.sites, n.trt = n.trt,
      n.days = n.days, n.years = n.years, sigma.sq = sigma.sq, lambda = lambda, ind = ind
    )
  )
}


