# This file provides a function for computing the design matrix for the 
# multivariate model described in Doser et al (2019). 

# Author: Jeffrey W. Doser


# computeDesignMatrix -----------------------------------------------------

computeDesignMatrix <- function(fileName, geoFileName, indexNames, yearsPostTrt) {
  
  # fileName = name of the csv file containing all data for the acoustic indices
  # geoFileName = name of the csv file containing the rain data (or other variable)
  # indexNames = character vector containing the names of acoustic indices used
  # yearsPostTrt = vector of the years in the model that are after the treatment
  #                occurred. 
  require(dplyr)
  require(tidyr)
  
  # Logit data transformation
  logit <- function(theta, a, b) {
    log((theta-a)/(b-theta))
  }
  
  # Inverse logit transformation
  logit.inv <- function(z, a, b) {
    b-(b-a)/(1+exp(z))
  }
  
  dat <- read.csv(fileName)
  rain <- read.csv(geoFileName, stringsAsFactors = FALSE)
  rain$siteInfo <- sub('.*/', '', rain$X)
  rain$year <- as.numeric(substr(rain$siteInfo, 5, 8))
  rain$day <- as.numeric(substr(rain$siteInfo, 11, 12))
  rain$time <- as.numeric(substr(rain$siteInfo, 14, 19))
  rain$X <- NULL
  rain$X.1 <- NULL
  rain$siteInfo <- NULL
  dat <- left_join(dat, rain, by = c("year", "day", "time"))
  
  # Transform variable to continuous for modeling with normal distribution
  for (i in 1:length(indexNames)) {
    if (indexNames[i] == 'ndsi') {
      dat <- dat %>% 
        mutate(ndsi = logit(ndsi, -1, 1))
    }
    if (indexNames[i] == 'H') {
      dat <- dat %>% 
        mutate(H = logit(H, 0, 1))
    } 
    if (indexNames[i] == 'aei') {
      dat <- dat %>% 
        mutate(aei = logit(aei, 0, 1))
    }
    if (indexNames[i] == 'aci') {
      dat <- dat %>% 
        mutate(aci = log(aci))
    }
    if (indexNames[i] == 'PSD1') {
      dat <- dat %>% 
        mutate(PSD1 = logit(PSD1, 0, 1))
    }
    if (indexNames[i] == 'PSD2') {
      dat <- dat %>% 
        mutate(PSD2 = logit(PSD2, 0, 1))
    }
    if (indexNames[i] == 'PSD3') {
      dat <- dat %>% 
        mutate(PSD3 = logit(PSD3, 0, 1))
    }
    if (indexNames[i] == 'PSD4') {
      dat <- dat %>% 
        mutate(PSD4 = logit(PSD4, 0, 1))
    }
    if (indexNames[i] == 'PSD5') {
      dat <- dat %>% 
        mutate(PSD5 = logit(PSD5, 0, 1))
    }
    if (indexNames[i] == 'PSD6') {
      dat <- dat %>% 
        mutate(PSD6 = logit(PSD6, 0, 1))
    }
    if (indexNames[i] == 'PSD7') {
      dat <- dat %>% 
        mutate(PSD7 = logit(PSD7, 0, 1))
    }
    if (indexNames[i] == 'PSD8') {
      dat <- dat %>% 
        mutate(PSD8 = logit(PSD8, 0, 1))
    }
    if (indexNames[i] == 'PSD9') {
      dat <- dat %>% 
        mutate(PSD9 = logit(PSD9, 0, 1))
    }
    if (indexNames[i] == 'PSD10') {
      dat <- dat %>% 
        mutate(PSD10 = logit(PSD10, 0, 1))
    }
  }
  
  # Transform data to long format
  long.dat <- gather(dat, key = "index", value = "value", indexNames)
  
  # Order data by each "individual" recording site/day/time combo
  ordered.long.dat <- long.dat %>% 
    arrange(recording.site, day, time, year)
  ordered.long.dat$X <- NULL
  n <- nrow(ordered.long.dat)
  
  # Create indicator variable for each "individual" soundscape
  ordered.long.dat$ind = 1
  for (i in 2:n) {
    if (sum(ordered.long.dat[i, c("month", "day", "time", "recording.site")] == 
            ordered.long.dat[i - 1, c("month", "day", "time", "recording.site")]) == 4) {
      ordered.long.dat$ind[i] = ordered.long.dat$ind[i-1]
    } else {
      ordered.long.dat$ind[i] = ordered.long.dat$ind[i-1] + 1
    }
  }
  
  
  # Construct design matrix
  n.indices <- length(indexNames)
  post.trt.years <- length(yearsPostTrt)
  X <- matrix(0, ncol = ((post.trt.years * 2) + 3) * n.indices, nrow = n)
  # X <- matrix(0, ncol = ((post.trt.years * 2) + 2) * n.indices, nrow = n)
  index.ind <- rep(1:n.indices, length.out = n)
  trt.ind <- ordered.long.dat$trt
  year <- ordered.long.dat$year
  num <- ncol(X) / n.indices
  for (i in 1:n.indices) {
    # Respective intercept for each psd value
    X[, (i-1)*num + 1] <- ifelse(index.ind == i, 1, 0)
    # Inherent difference parameter for each index
    X[, (i-1)*num + 2] <- ifelse(index.ind == i & trt.ind == 1, 1, 0)
    # Year effects
    X[, (i-1)*num + 3] <- ifelse(index.ind == i & year == yearsPostTrt[1], 1, 0)
    X[, (i-1)*num + 4] <- ifelse(index.ind == i & year == yearsPostTrt[2], 1, 0)
    X[, (i-1)*num + 5] <- ifelse(index.ind == i & year == yearsPostTrt[3], 1, 0)
    X[, (i-1)*num + 6] <- ifelse(index.ind == i & year == yearsPostTrt[4], 1, 0)
    X[, (i-1)*num + 7] <- ifelse(index.ind == i & year == yearsPostTrt[5], 1, 0)
    # Treatment effects for each year
    X[, (i-1)*num + 8] <- ifelse(index.ind == i & year == yearsPostTrt[1] & trt.ind == 1, 1, 0)
    X[, (i-1)*num + 9] <- ifelse(index.ind == i & year == yearsPostTrt[2] & trt.ind == 1, 1, 0)
    X[, (i-1)*num + 10] <- ifelse(index.ind == i & year == yearsPostTrt[3] & trt.ind == 1, 1, 0)
    X[, (i-1)*num + 11] <- ifelse(index.ind == i & year == yearsPostTrt[4] & trt.ind == 1, 1, 0)
    X[, (i-1)*num + 12] <- ifelse(index.ind == i & year == yearsPostTrt[5] & trt.ind == 1, 1, 0)
    # Rain effect
    X[, (i-1)*num + 13] <- ifelse(index.ind == i, ordered.long.dat$rain, 0)
  }
  
  return(list(
    X = X, dat = ordered.long.dat, n.indices = n.indices
  ))
}
