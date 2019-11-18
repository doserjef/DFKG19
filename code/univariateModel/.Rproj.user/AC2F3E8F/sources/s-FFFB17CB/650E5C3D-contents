rm(list = ls())
library(coda)
library(ggplot2)
library(tidyr)
library(dplyr)

# Posterior summary information from the single index model output, along
# with some posterior visualizations of important model results

# Author: Jeffrey W. Doser
# summary.R ---------------------------------------------------------------

n.iter <- 500
theta.1 <- t(matrix(scan("theta-samples-1"), ncol = n.iter, byrow = F))

n.alpha <- 13
# n.alpha <- 12
# n.alpha <- 11
colnames(theta.1) <- c(paste0("alpha.", 0:(n.alpha-1)), "tau.sq.", "sigma.sq")


burn.in <- floor(0.1 * n.iter)
sub <- (burn.in+1):n.iter

# Post burn-in quantiles for each parameter
m.quants.1 <- summary(window(mcmc(theta.1), start = burn.in))$quantiles
m.quants.1

# Assess convergence with multiple chains if desired ----------------------
# Would first need to run the model with different starting values and save
# with different file names as indicated by the commented out code below
# theta.2 <- t(matrix(scan("theta-samples-2"), ncol = n.iter, byrow =F))
# theta.3 <- t(matrix(scan("theta-samples-3"), ncol = n.iter, byrow = F))
# colnames(theta.2) <- colnames(theta.1)
# colnames(theta.3) <- colnames(theta.2)
# m.quants.2 <- summary(window(mcmc(theta.2), start = burn.in))$quantiles
# m.quants.3 <- summary(window(mcmc(theta.3), start = burn.in))$quantiles
# chains <- mcmc.list(mcmc(theta.1), mcmc(theta.2), mcmc(theta.3))
# plot(chains, density = FALSE)
# gelman.diag(chains)

# if overlaps with zero, suggests there are no differences in b/w control and treatment. 
ctrl.alpha.post <- as.data.frame(m.quants.1[3:7, ])
row.names(ctrl.alpha.post) <- c("ctrl2014", "ctrl2015", "ctrl2016", "ctrl2017", "ctrl2018")
years <- c(2014:2018)
ctrl.alpha.post$years <- years
trt.alpha.post <- as.data.frame(m.quants.1[8:12, ])
row.names(trt.alpha.post) <- c("trt2014", "trt2015", "trt2016", "trt2017", "trt2018")
trt.alpha.post$years <- years

single.plot.dat <- data.frame(year = trt.alpha.post$years, 
                              low = trt.alpha.post[, 1], 
                              med = trt.alpha.post[, 3], 
                              high = trt.alpha.post$`97.5%`, 
                              index = 'NDSI')

# Figure to summarize treatment effects -----------------------------------

png("../../figures/ndsiYearEffects.png", height = 480, width = 580)
theme_set(theme_bw(base_size = 18))
ggplot(data = single.plot.dat, aes(x = year)) + 
  geom_point(aes(y = med), size = 2.5) + 
  geom_segment(aes(x = year, y = low, xend = year, yend = high), 
               color = 'red', size = 0.9) + 
  facet_wrap(index ~ ., scales = "free_y") + 
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'black') + 
  labs(x = 'Year', y = 'Year Effect') + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18),
        plot.margin = margin(t = 16, r = 16, b = 0, l = 5, unit = "pt"))
dev.off()


# Facet plot of Gibbs sampler output --------------------------------------

vals.post.burn <- window(mcmc(theta.1), start = burn.in)

dat <- as.data.frame(vals.post.burn)

dat.long <- gather(dat, key = 'Variable', value = 'val')

dat.long.sub <- dat.long %>% 
  filter(Variable %in% c("alpha.1", "alpha.7", "alpha.8", "alpha.9", "alpha.10", 
                         "alpha.11", "alpha.12", "sigma.sq", "tau.sq."))

dat.long.sub$Variable <- factor(dat.long.sub$Variable, 
                                levels = c('alpha.1', 'alpha.7', 'alpha.8', 
                                           'alpha.9', 'alpha.10', 'alpha.11', 
                                           'alpha.12', 'sigma.sq', 'tau.sq.'))

# Code up names to enable inclusion of math symbols in the facet plot titles
plot_names <- c(
  `alpha.1` = 'scriptstyle(bgroup("", a, ")"))~alpha[2]',
  `alpha.7` = 'scriptstyle(bgroup("", b, ")"))~alpha[8]',
  `alpha.8` = 'scriptstyle(bgroup("", c, ")"))~alpha[9]',
  `alpha.9` = 'scriptstyle(bgroup("", d, ")"))~alpha[10]', 
  `alpha.10` = 'scriptstyle(bgroup("", e, ")"))~alpha[11]', 
  `alpha.11` = 'scriptstyle(bgroup("", f, ")"))~alpha[12]', 
  `alpha.12` = 'scriptstyle(bgroup("", g, ")"))~alpha[13]', 
  `sigma.sq` = 'scriptstyle(bgroup("", h, ")"))~sigma^2', 
  `tau.sq.` = 'scriptstyle(bgroup("", i, ")"))~tau^2'
)

png("../../figures/univariateGibbsOut.png", width = 960, height = 480)
theme_set(theme_bw(base_size = 18))
ggplot(data = dat.long.sub, aes(x = val, fill = 'black')) + 
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = 'dashed', col = 'black') +
  guides(fill = FALSE) + 
  facet_wrap(Variable ~ ., scales = "free_y", labeller = as_labeller(plot_names, label_parsed)) + 
  labs(x = 'Parameter Value', y = 'Density') + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18), 
        plot.margin = margin(t = 16, r = 16, b = 0, l = 5, unit = "pt")) 
dev.off()
