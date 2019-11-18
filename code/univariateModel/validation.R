rm(list = ls())
library(coda)
library(ggplot2)
library(tidyr)
library(dplyr)

# Code to compute and display model validation results from Doser et al (2019): 
# Assessing soundscape disturbance through hierarchical modelsand acoustic indices: 
# a case study on a shelterwood logged northern Michigan forest

# Author: Jeffrey W. Doser
# Validation  -------------------------------------------------------------

n.iter <- 500
theta.full <- t(matrix(scan("theta-samples-1"), ncol = n.iter, byrow = F))
n.alpha <- 13
colnames(theta.full) <- c(paste0("alpha.", 0:(n.alpha - 1)), "tau.sq.", "sigma.sq")
theta.no.diffs <- t(matrix(scan("theta-no-diffs-samples"), ncol = n.iter, byrow = F))
colnames(theta.no.diffs) <- colnames(theta.full)[-2]
theta.no.rain <- t(matrix(scan("theta-no-rain-samples"), ncol = n.iter, byrow = F))
colnames(theta.no.rain) <- colnames(theta.full)[-13]
theta.no.re <- t(matrix(scan("no-beta-theta-samples"), ncol = n.iter, byrow = F))
colnames(theta.no.re) <- colnames(theta.full)[-14]
theta.no.all <- t(matrix(scan("no-beta-diff-rain-theta-samples"), ncol = n.iter, 
                         byrow = F))
colnames(theta.no.all) <- colnames(theta.full)[-c(2, 13, 14)]

burn.in <- floor(0.1 * n.iter)
sub <- (burn.in+1):n.iter

quants.full <- summary(window(mcmc(theta.full), start = burn.in))$quantiles
quants.no.diffs <- summary(window(mcmc(theta.no.diffs), start = burn.in))$quantiles
quants.no.rain <- summary(window(mcmc(theta.no.rain), start = burn.in))$quantiles
quants.no.re <- summary(window(mcmc(theta.no.re), start = burn.in))$quantiles
quants.no.all <- summary(window(mcmc(theta.no.all), start = burn.in))$quantiles

all.quants <- rbind(quants.full, quants.no.diffs, quants.no.rain, 
                    quants.no.re, quants.no.all)
vals <- row.names(all.quants)
all.quants <- as.data.frame(all.quants, row.names = FALSE)
all.quants$vals <- vals

years <- c(2014:2018)
trt.alpha.post <- all.quants %>% 
  filter(vals %in% c('alpha.7', 'alpha.8', 'alpha.9', 'alpha.10', 'alpha.11'))
# trt.alpha.post$model <- rep(c('Full', paste('No ', expression(alpha[2]), sep = ''),
#                               paste('No ', expression(alpha[13]), sep = ''),
#                               paste('No ', expression(tau^2), sep = ''),
#                               paste('No ', expression(alpha[2]), ', ', 
#                                     expression(alpha[13]), ', ', 
#                                     expression(tau^2), sep = '')), 
#                             each = length(years))
trt.alpha.post$model <- factor(rep(c('Full', 'No Inherent Differences', 'No Rain Effect', 
                                     'No Random Effects', 'Basic'), 
                               each = length(years)))
trt.alpha.post$model <- factor(rep(1:5, each = length(years)), levels = c(1, 2, 3, 4, 5), 
                               labels = c('Full', 'No Inherent Differences', 
                                          'No Rain Effect', 'No Random Effects', 'Basic'))
trt.alpha.post$years <- rep(years, 5)
trt.alpha.post$years[1:5] <- trt.alpha.post$years[1:5] - .2
trt.alpha.post$years[6:10] <- trt.alpha.post$years[6:10] - .1
trt.alpha.post$years[11:15] <- trt.alpha.post$years[11:15]
trt.alpha.post$years[16:20] <- trt.alpha.post$years[16:20] + .1
trt.alpha.post$years[21:25] <- trt.alpha.post$years[21:25] + .2

single.plot.dat <- data.frame(year = trt.alpha.post$years, 
                              low = trt.alpha.post[, 1], 
                              med = trt.alpha.post[, 3], 
                              high = trt.alpha.post$`97.5%`, 
                              vals = trt.alpha.post$vals,
                              model = trt.alpha.post$model,
                              index = 'NDSI')

png("../../figures/ndsiValidation.png", height = 480, width = 780)
theme_set(theme_bw(base_size = 18))
ggplot(data = single.plot.dat, aes(x = year, col = model, shape = model)) + 
  geom_point(aes(y = med), size = 2.5) + 
  geom_segment(aes(x = year, y = low, xend = year, yend = high), size = 0.9) + 
  facet_wrap(index ~ ., scales = "free_y") + 
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'black') + 
  labs(x = 'Year', y = 'Year Effect', col = 'Model', shape = 'Model') + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18),
        plot.margin = margin(t = 16, r = 16, b = 0, l = 5, unit = "pt"))
dev.off()

