rm(list = ls())
library(dplyr)

# This document creates the design matrix for a univariate model. 
# Here the NDSI is used as the single acoustic index.

# Author: Jeffrey W. Doser
# 03-yearDesignMatrix -----------------------------------------------------

# Read in psd/ndsi data
dat <- read.csv("../data/full-ndsi-psd-ordered.csv")
# Read in rain data that was manually determined from the recordings
rain <- read.csv("../data/geophony.csv", stringsAsFactors = FALSE)
# Add information to the rain data frame from the recording name
rain$siteInfo <- sub('.*/', '', rain$X)
rain$year <- as.numeric(substr(rain$siteInfo, 5, 8))
rain$day <- as.numeric(substr(rain$siteInfo, 11, 12))
rain$time <- as.numeric(substr(rain$siteInfo, 14, 19))
rain$X <- NULL
rain$X.1 <- NULL
rain$siteInfo <- NULL

# Join the acoustic index data with the rain data
dat <- left_join(dat, rain, by = c("year", "day", "time"))

n <- nrow(dat)
# Each unique site/day/time combination is an individual
# Create individual identity (i.e., assign each individual soundscape a number)
n.ind.counter <- 1
dat$ind = 1
for (i in 2:n) {
  if (sum(dat[i, c("month", "day", "time", "recording.site")] == 
          dat[i - 1, c("month", "day", "time", "recording.site")]) == 4) {
    dat$ind[i] = dat$ind[i-1]
  } else {
    dat$ind[i] = dat$ind[i-1] + 1
  }
}

# Number of recording years after the treatment
post.trt.years <- 5
X <- matrix(0, nrow = n, ncol = 3 + post.trt.years * 2)

# Create design matrix ----------------------------------------------------
# General intercept
X[, 1] <- 1
# Intercept for treatment sites only. Represents inherent differences b/w 
# control and treatment sites. Treatment sites given value 1, control sites
# given value 0
X[, 2][dat$trt == 1] <- 1
# Overall year effects regardless of site. Take value 1 if in the respective 
# year, value 0 if not in that year
X[, 3] <- ifelse(dat$year == 2014, 1, 0)
X[, 4] <- ifelse(dat$year == 2015, 1, 0)
X[, 5] <- ifelse(dat$year == 2016, 1, 0)
X[, 6] <- ifelse(dat$year == 2017, 1, 0)
X[, 7] <- ifelse(dat$year == 2018, 1, 0)
# Treatment effect sites in each year. Takes value 1 if in the respective
# year and a treatment site, 0 if otherwise
X[, 8] <- ifelse(dat$year == 2014 & dat$trt == 1, 1, 0)
X[, 9] <- ifelse(dat$year == 2015 & dat$trt == 1, 1, 0)
X[, 10] <- ifelse(dat$year == 2016 & dat$trt == 1, 1, 0)
X[, 11] <- ifelse(dat$year == 2017 & dat$trt == 1, 1, 0)
X[, 12] <- ifelse(dat$year == 2018 & dat$trt == 1, 1, 0)
# Influence of rain. 1 if there is rain, 0 if not
X[, 13] <- dat$rain

# Read out the design matrix
write.table(X, "../data/X-year", row.names = FALSE, col.names = FALSE, sep = "\t")
# Read out the ndsi values
write.table(dat$ndsi, "../data/ndsi", row.names = FALSE, col.names = FALSE, sep = "\t")
# Read out the individual soundscape identifier
write.table(dat$ind, "../data/ind", row.names = FALSE, col.names = FALSE, sep = "\t")
