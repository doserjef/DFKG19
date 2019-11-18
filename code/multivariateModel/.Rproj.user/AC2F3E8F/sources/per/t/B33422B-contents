rm(list = ls())
source("sim.R")
source("sim-multivariateSoundModel.R")
library(coda)

# This is the main program for testing the multivariate model with simulated
# data

# Author: Jeffrey W. Doser

# Main Program ------------------------------------------------------------

# Simulate the data
dat <- sim.multi.dat()
n.iter <- 300
out <- sim.multivar.sound.model(n.iter = n.iter, 
                                X = dat$X, 
                                y = dat$y, 
                                n.sites = dat$n.sites, 
                                n.indices = nrow(dat$lambda), 
                                ind = dat$ind, 
                                alpha.start = 0, 
                                beta.start = 0, 
                                sigma.sq.start = 1,
                                lambda.start = 1)

sigma.sq.samples <- out$sigma.sq.samples
lambda.samples <- out$lambda.samples
alpha.samples <- out$alpha.samples
beta.samples <- out$beta.samples
# Summary -----------------------------------------------------------------

burn.in <- floor(0.5 * n.iter)
sub <- (burn.in+1):n.iter

# Look at process varaince
plot(window(mcmc(sigma.sq.samples), start = burn.in), density = FALSE)
summary(window(mcmc(sigma.sq.samples), start = burn.in))
dat$sigma.sq

# Lambda covariance matrix
summary(window(mcmc(t(lambda.samples)), start = burn.in))

# Burn in. Can change the 0.5 to change the amount of burn.in required
burn.in <- floor(0.5 * n.iter)
sub <- (burn.in+1):n.iter

# Fixed coefficients
alpha.hat.mean <- apply(alpha.samples[, sub], 1, mean)
alpha.true <- dat$alpha
plot(alpha.true, alpha.hat.mean, pch = 19)
lines(alpha.true, alpha.true, col = 'red')
summary(mcmc(window(t(alpha.samples), start = burn.in)))

summary(mcmc(window(t(lambda.samples), start = burn.in)))


# Individual random effects
beta.hat.mean <- apply(beta.samples[, sub], 1, mean)
beta.true <- dat$beta
plot(beta.true, beta.hat.mean, pch = 19)
lines(beta.true, beta.true, col = 'red')
