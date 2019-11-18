rm(list = ls())
source("sim.R")
source("singleIndexModel.R")
library(coda)

# Main program for testing the univariate model with simulated data

# Author: Jeffrey W. Doser

# sim-main.R --------------------------------------------------------------

# Simulate data
dat <- sim.univariate.data()
# Run the model
out <- univar.sound.model(n.iter = 500, n.sites = dat$n.sites, data = dat$y, 
                          ind = dat$ind, X = dat$X, alpha.start = 0, 
                          beta.start = 0, sigma.sq.start = 1, 
                          tau.sq.start = 1,index = 'sim')

alpha.samples <- out$alpha.samples
tau.sq.samples <- out$tau.sq.samples
sigma.sq.samples <- out$sigma.sq.samples
beta.samples <- out$beta.samples
n.iter <- out$n.iter

# Summary -----------------------------------------------------------------

plot(mcmc(t(alpha.samples)),density = FALSE)
plot(mcmc(tau.sq.samples), density = FALSE)
plot(mcmc(sigma.sq.samples), density = FALSE)

theta <- cbind(t(alpha.samples), tau.sq.samples, sigma.sq.samples)
n.alpha <- nrow(alpha.samples)
colnames(theta) <- c(paste0("alpha.", 1:n.alpha), "tau.sq.", "sigma.sq")

burn.in <- floor(0.5 * n.iter)
sub <- (burn.in+1):n.iter

# Determine if true simulated values are inside the 95% credible interval for
# the estimated parameters
summary(window(mcmc(theta), start = burn.in))
dat$alpha
dat$sigma.sq
dat$tau.sq
