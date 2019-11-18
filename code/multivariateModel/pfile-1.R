alpha <- rep(0, n.alpha)
beta <- matrix(0, nrow = n.ind * n.psd, ncol = n.beta)
sigma.sq <- 1
lambda <- diag(.1, n.psd)

num <- 1