rm(list = ls())
source("multivarSoundModel.R")
library(coda)
library(corrplot)
library(ggplot2)
library(dplyr)

# This is the main file for running the multivariate model described in 
# Doser et al (2019) Assessing soundscape disturbance through hierarchical models
# and acoustic indices: a case study on a shelterwood logged northern Michigan forest. 
# Also includes code for computation of figures for analysis of model. 

# Note code for convergence assessment is not included in this file, but can be 
# implemented in a similar manner as displayed in the univariateModel using 
# the Gelman-Rubin diagnostic. 

# Author: Jeffrey W. Doser

# Run model ---------------------------------------------------------------
n.iter <- 500
# Change to adhere to your directory structure. 
fileName <- "../../data/orderedMultivariateData.csv"
geoFileName <- "../../data/geophony.csv"
indexNames <- c('H', 'aci', 'ndsi', 'aei', paste('PSD', 1:10, sep = ''))
yearsPostTrt <- 2014:2018
alpha.start <- 1
beta.start <- 0
sigma.sq.start <- 1
lambda.start <- 1
out1 <- multivar.sound.model(n.iter = n.iter, 
                             fileName,
                             geoFileName, 
                             indexNames,
                             yearsPostTrt, 
                             alpha.start, 
                             beta.start, 
                             sigma.sq.start, 
                             lambda.start)

alpha.samples <- out1$alpha.samples
beta.samples <- out1$beta.samples
sigma.sq.samples <- out1$sigma.sq.samples
lambda.samples <- out1$lambda.samples


# Summary -----------------------------------------------------------------

burn.in <- floor(0.1 * n.iter)
sub <- (burn.in+1):n.iter

# alpha 
alpha.vals <- summary(window(mcmc(t(alpha.samples)), start = burn.in))$quantiles

# sigma.sq
summary(window(mcmc(sigma.sq.samples)), start = burn.in)

#lambda
vals <- summary(window(mcmc(t(lambda.samples)), start = burn.in))$quantiles
n.indices <- length(indexNames)
lambda.med.mat <- matrix(vals[, 3], nrow = n.indices)
cov2cor(lambda.med.mat)
lambda.lower.mat <- matrix(vals[, 1], nrow = n.indices)
lambda.upper.mat <- matrix(vals[, 5], nrow = n.indices)
cov2cor(lambda.lower.mat)
cov2cor(lambda.upper.mat)

# Visualize the random effects correlation matrix -------------------------
corr.med <- cov2cor(lambda.med.mat)
colnames(corr.med) <- toupper(indexNames)
rownames(corr.med) <- toupper(indexNames)
corr.low <- cov2cor(lambda.lower.mat)
colnames(corr.low) <- toupper(indexNames)
rownames(corr.low) <- toupper(indexNames)
corr.high <- cov2cor(lambda.upper.mat)
colnames(corr.high) <- toupper(indexNames)
rownames(corr.high) <- toupper(indexNames)

sig.mat <- matrix(0, nrow(corr.med), ncol(corr.med))
# Note that 0 denotes "significant", 1 denotes "insignificant
sig.mat <- ifelse((0 > corr.low) & (0 < corr.high), 1, 0)


col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
png("../../figures/correlationDiagram.png")
corrplot(corr.med, method="color", col=col(200),  
         type="upper", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # hide correlation coefficient on the principal diagonal
         diag=FALSE, 
         # "Significance"
         p.mat = sig.mat, sig.level = 0.5, insig = 'blank'
)
dev.off()

# Rain Effect
alpha.out.sub <- alpha.vals[13:182, ]
rain.out <- alpha.out.sub[seq(1, nrow(alpha.out.sub), 13), ]

#  Check out the non-treatment effects ------------------------------------
alpha.nontrt.out <- alpha.vals[2:171, ]
nontrt.out <- alpha.nontrt.out[seq(1, nrow(alpha.nontrt.out), 13), ]


# Year Effects for Each Acoustic Index in a Facet Plot --------------------

# A single plot for a single index
rownames(alpha.vals) <- paste(rep(indexNames, each = 13), 1:13, sep = "")
plot.dat <- data.frame(low = alpha.vals[, 1], 
                       med = alpha.vals[, 3], 
                       high = alpha.vals[, 5], 
                       alpha = rep(1:13, n.indices), 
                       index = rep(indexNames, each = 13))

plot.trt.dat <- plot.dat %>% 
  filter(alpha %in% c(8:12))
plot.trt.dat$year <- rep(2014:2018, n.indices)

plot.trt.dat$index <- toupper(plot.trt.dat$index)
plot.trt.dat$index <- factor(plot.trt.dat$index, 
                             levels = c('H', 'ACI', 'NDSI', 'AEI', 
                                        "PSD1", "PSD2", "PSD3", 'PSD4', 
                                        'PSD5', 'PSD6', 'PSD7', 'PSD8', 
                                        'PSD9', 'PSD10'))

plot.names <- c(
  `H` = 'scriptstyle(bgroup("", a, ")"))~H', 
  `ACI` = 'scriptstyle(bgroup("", b, ")"))~ACI',
  `NDSI` = 'scriptstyle(bgroup("", c, ")"))~NDSI',
  `AEI` = 'scriptstyle(bgroup("", d, ")"))~AEI', 
  `PSD1` = 'scriptstyle(bgroup("", e, ")"))~PSD[1]',
  `PSD2` = 'scriptstyle(bgroup("", f, ")"))~PSD[2]',
  `PSD3` = 'scriptstyle(bgroup("", g, ")"))~PSD[3]',
  `PSD4` = 'scriptstyle(bgroup("", h, ")"))~PSD[4]',
  `PSD5` = 'scriptstyle(bgroup("", i, ")"))~PSD[5]',
  `PSD6` = 'scriptstyle(bgroup("", j, ")"))~PSD[6]',
  `PSD7` = 'scriptstyle(bgroup("", k, ")"))~PSD[7]',
  `PSD8` = 'scriptstyle(bgroup("", l, ")"))~PSD[8]',
  `PSD9` = 'scriptstyle(bgroup("", m, ")"))~PSD[9]',
  `PSD10` = 'scriptstyle(bgroup("", n, ")"))~PSD[10]'
)

png("../../figures/multivariateTreatmentEffects.png", height = 780, width = 650)
theme_set(theme_bw(base_size = 18))
ggplot(data = plot.trt.dat, aes(x = year)) + 
  geom_point(aes(y = med)) + 
  geom_segment(aes(x = year, y = low, xend = year, yend = high), 
               color = 'red') + 
  facet_wrap(index ~ ., scales = "free_y", ncol = 3, 
             labeller = as_labeller(plot.names, label_parsed)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'black') + 
  labs(x = 'Year', y = 'Year Effect') + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18),
        plot.margin = margin(t = 16, r = 16, b = 0, l = 5, unit = "pt"))
dev.off()
