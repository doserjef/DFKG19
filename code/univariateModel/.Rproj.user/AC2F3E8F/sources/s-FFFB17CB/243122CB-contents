rm(list = ls())
library(coda)

# This program runs the model without any random effects. Used in the 
# model comparison using WAIC portion of Doser et al (2019).

# Author: Jeffrey W. Doser

# Main Program ------------------------------------------------------------

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

n.samples <- 500
n.sites <- 13
# ndsi <- as.matrix(read.table("../../data/ndsi"))[, 1]
ndsi <- as.matrix(read.csv("../../data/orderedIndices.csv")$ndsi)
y <- as.matrix(logit(ndsi, a = -1, b = 1))
X <- as.matrix(read.table("../../data/X-year"))

# For model validation, can remove different important parameters and see 
# the effects this has on inference.

# Remove the inherent differences parameter
# X <- X[, -2]

# Remove the rain parameter
# X <- X[, -13]

# Remove rain and inherent differences
# X <- X[, -c(2, 13)]

ind <- read.table("../../data/ind")[, 1]
n.ind <- length(unique(ind))
n <- nrow(X)
n.alpha <- dim(X)[2]
s <- rep(0, n.ind)
for (i in 1:n.ind) {
  curr <- sum(ind == i)
  s[i] <- curr
}

# Priors ------------------------------------------------------------------

alpha.mu <- rep(0, n.alpha)
alpha.var <- rep(10000, n.alpha)

sigma.sq.a <- 2
sigma.sq.b <- 1



# Starting values ---------------------------------------------------------

alpha <- rep(0, n.alpha)
sigma.sq <- 1
y.hat <- rep(0, n)


# Sampler Prep ------------------------------------------------------------

alpha.samples <- matrix(0, nrow = n.alpha, ncol = n.samples)
fitted.samples <- matrix(0, nrow = n, ncol = n.samples)
sigma.sq.samples <- rep(0, n.samples)
post.dens.samples <- matrix(0, n, n.samples)
post.log.dens.samples <- matrix(0, n, n.samples)
XTX <- t(X)%*%X

# Gibbs Sampler -----------------------------------------------------------

for (k in 1:n.samples) {
  
  # Update alpha samples
  V <- chol2inv(chol(XTX / sigma.sq + diag(1/alpha.var, n.alpha)))
  sum.diffs <- rep(0, n.alpha)
  for (i in 1:n.ind) {
    X.i <- X[ind == i, ]
    y.i <- y[ind == i]
    # W.i <- W[ind == i, ]
    sum.diffs <- sum.diffs + t(X.i) %*% (y.i)
  }
  v <- sum.diffs / sigma.sq
  alpha <- rmvn(1, V%*%v, V)
  
  # Update sigma.sq samples
  a <- sigma.sq.a + .5 * sum(s)
  b.sum <- 0
  for (i in 1:n.ind) {
    X.i <- X[ind == i, ]
    y.i <- y[ind == i]
    # W.i <- W[ind == i, ]
    b.sum <- b.sum + t((y.i - (X.i%*%alpha))) %*% 
      (y.i - (X.i%*%alpha))
  }
  b <- sigma.sq.b + .5 * b.sum
  sigma.sq <- rigamma(1, a, b)
  
  # Fitted values
  curr.length <- 0 
  for (i in 1:n.ind) {
    x.i <- X[ind == i, ]
    # W.i <- W[ind == i, ]
    V <- diag(sigma.sq, s[i])
    mu <- x.i %*% alpha
    y.hat[(curr.length + 1):(curr.length + s[i])] <- rmvn(n = 1, mu = mu, V = V)
    curr.length <- curr.length + s[i]
  }
  
  # Get values for WAIC ---------------------------------------------------
  mu <- X %*% alpha
  post.dens.samples[, k] <- dnorm(y, mu, sigma.sq)
  post.log.dens.samples[, k] <- dnorm(y, mu, sigma.sq, log = TRUE)
  
  # Save Samples
  alpha.samples[, k] <- alpha
  sigma.sq.samples[k] <- sigma.sq
  fitted.samples[, k] <- y.hat
  
  print(paste(k/n.samples * 100, " percent complete", sep = ""))
  
}

# Compute WAIC
ppd <- sum(log(apply(post.dens.samples, 1, sum) / n.samples))
p.d <- sum(apply(post.log.dens.samples, 1, var))
waic <- -2 * ppd + 2 * p.d


# Summary -----------------------------------------------------------------

theta <- cbind(t(alpha.samples), sigma.sq.samples)
colnames(theta) <- c(paste0("alpha.", 0:(n.alpha-1)), "sigma.sq")


write.table(theta, "no-beta-diff-rain-theta-samples", col.names = FALSE, 
            row.names = FALSE, sep = "\t")
write.table(fitted.samples, "no-beta-diff-rain-fitted-samples", col.names = FALSE, 
            row.names = FALSE, sep = "\t")


