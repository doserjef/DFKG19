sim.univariate.data <- function(n.sites = 13, n.dc = 4, n.trt = 4,
                                n.days = 30, n.years = 5, 
                                alpha = c(1, -1, 0.2, 0.1, -.1, -.8, -.5, -.6, -.3, 1), 
                                sigma.sq = 2, tau.sq = 0.3, n.beta = 1) {
  
  # n.sites = number of recording sites
  # n.dc = number of recordings on each recording day
  # n.trt = number of treatment sites (< n.sites)
  # n.days = number of recording days
  # n.years = number of recording years
  # alpha = fixed effects regression coefficients
  # sigma.sq = process error
  # tau.sq = random effects variance
  # n.beta = set to 1 for individual random effects, if > 1 means you're adding in 
  #          additional types of random effects
  
  # This function simulates data 
  # Author: Jeffrey W. Doser
  
  # Multivariate normal random number generator -----------------------------
  
  rmvn <- function(n, mu=0, V = matrix(1)){
    p <- length(mu)
    if(any(is.na(match(dim(V),p))))
      stop("Dimension problem!")
    D <- chol(V)
    t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
  }
  
  # Logit data transformation
  logit <- function(theta, a, b) {
    log((theta-a)/(b-theta))
  }
  
  # Inverse logit transformation
  logit.inv <- function(z, a, b) {
    b-(b-a)/(1+exp(z))
  }
  
  
  # Simulate Data -----------------------------------------------------------
  
  set.seed(1)
  n <- n.sites * n.dc * n.days * n.years
  n.ind <- n.days * n.dc * n.sites
  
  beta <- matrix(rnorm(n.ind * n.beta, 0, sd = sqrt(tau.sq)), nrow = n.ind, ncol = n.beta)
  # beta.1 <- rnorm(n.ind, 0, sd = sqrt(tau.1))
  # beta <- cbind(beta.0, beta.1)
  
  n.ctrl <- n.sites - n.trt
  
  
  # Construct design matrix -------------------------------------------------
  
  # Create an indicator for treatment groups
  # Note that treatment groups are listed first
  trt.ind <- rep(0, n)
  # Note this assumes the treatment occurred in the second year
  # i.e., recordings are obtained one year before the treatment and continue 
  # for n.years - 1 more years.
  trt.ind[1:(600*4)] <- 1
  year <- rep(c(2013, 2014, 2015, 2016, 2017), length.out = n)
  # Intercept
  x.0 <- rep(1, n) 
  # Difference b/w treatment and control
  x.1 <- rep(0, n)
  x.1[trt.ind == 1] <- 1
  # # Treatment effect on NDSI
  # x.2 <- rep(0, n)
  # x.2[trt.ind == 1 & year > 2013] <- 1
  # Year effects for all sites
  x.2 <- ifelse(year == 2014, 1, 0)
  x.3 <- ifelse(year == 2015, 1, 0)
  x.4 <- ifelse(year == 2016, 1, 0)
  x.5 <- ifelse(year == 2017, 1, 0)
  # Year effects of treatment
  x.6 <- ifelse(year == 2014 & trt.ind == 1, 1, 0)
  x.7 <- ifelse(year == 2015 & trt.ind == 1, 1, 0)
  x.8 <- ifelse(year == 2016 & trt.ind == 1, 1, 0)
  x.9 <- ifelse(year == 2017 & trt.ind == 1, 1, 0)
  # Construct matrix
  X <- cbind(x.0, x.1, x.2, x.3, x.4, x.5, x.6, x.7, x.8, x.9)
  
  # Indicator for individual
  ind <- rep(1:n.ind, each = n.years)
  
  # Random effects design matrix
  # W <- X[, 1:2]
  W <- as.matrix(X[, 1])
  
  # Construct likelihood
  V <- diag(sigma.sq, n.years)
  y <- rep(0, n)
  for (i in 1:n.ind) {
    x.i <- X[ind == i, ]
    # W.i <- W[ind == i, ]
    W.i <- as.matrix(W[ind == i, ])
    # v <- x.i %*% alpha + rep(1, n.years)*beta.0[i] + w.i * beta.1[i]
    v <- x.i %*% alpha + W.i%*%beta[i, ]
    y[((i-1)*n.years + 1):(((i-1)*n.years)+n.years)] <- rmvn(n = 1, mu = v, V = V)
  }
  
  return(
    list(
      X = X, y = y, W = W, beta = beta, n.sites = n.sites,
      n.trt = n.trt, n.days = n.days, n.years = n.years, alpha = alpha, 
      sigma.sq = sigma.sq, tau.sq = tau.sq, n.beta = n.beta, ind = ind
    )
  )
}