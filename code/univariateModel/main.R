rm(list = ls())
library(coda)
source("singleIndexModel.R")

# Main program to call univariate acoustic index model and output the 
# posterior samples for all parameters of interest.

# Author: Jeffrey W. Doser

# Main Program ------------------------------------------------------------


out <- univar.sound.model(n.iter = 500, 
                          n.sites = 13, 
                          data = read.csv("../../data/orderedIndices.csv")$ndsi, 
                          ind = read.table("../../data/ind")[, 1], 
                          X = as.matrix(read.table("../../data/X-year")), 
                          alpha.start = 0,
                          beta.start = 0,
                          sigma.sq.start = 1,
                          tau.sq.start = 1, 
                          index = 'NDSI')

alpha.samples <- out$alpha.samples
sigma.sq.samples <- out$sigma.sq.samples
tau.sq.samples <- out$tau.sq.samples
post.dens.samples <- out$post.dens.samples
post.log.dens.samples <- out$post.log.dens.samples
n.iter <- out$n.iter

# Compute WAIC for model comparison
ppd <- sum(log(apply(post.dens.samples, 1, sum) / n.iter))
p.d <- sum(apply(post.log.dens.samples, 1, var))
waic <- -2 * ppd + 2 * p.d


# Summary -----------------------------------------------------------------

theta <- cbind(t(alpha.samples), tau.sq.samples, sigma.sq.samples)
n.alpha <- nrow(alpha.samples)
colnames(theta) <- c(paste0("alpha.", 0:(n.alpha-1)), "tau.sq.", "sigma.sq")
num <- 1

# Output all posterior samples for analysis in summary.R
write.table(theta, paste("theta-samples", num, sep = "-"), col.names = FALSE, 
            row.names = FALSE, sep = "\t")


