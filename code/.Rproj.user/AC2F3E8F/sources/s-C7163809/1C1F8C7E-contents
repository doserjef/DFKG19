rm(list = ls())
source("soundscape.R")
library(dplyr)

# Compute additional acoustic indices used in the multivariate model
# in Doser et al (2019) preprint. 
# 08-additionalIndexComputations ------------------------------------------

# Compute acoustic indices one at a time using multiple_sounds function. 
# Do this for each index of interest.
# multiple_sounds(directory = "../recordings/", resultfile = "../data/aci.csv",
                # soundindex ="acoustic_complexity", no_cores = 2)


# Organize Data sets ------------------------------------------------------

# H
full.dat <- extractMultipleSounds("../data/H.csv", index = "H")

# Properly structure data
n <- dim(full.dat)[1]
ordered.dat <- full.dat %>% 
  arrange(recording.site, day, time, year)
head(ordered.dat, 30)

# ACI
aci.dat <- extractMultipleSounds("../data/aci.csv", index = 'ACI')

aci.ordered <- aci.dat %>% 
  arrange(recording.site, day, time, year) %>% 
  select(ACI)

ordered.dat$aci <- aci.ordered$ACI

# NDSI
ndsi.dat <- extractMultipleSounds("../data/ndsi.csv", index = 'ndsi')
ndsi.ordered <- ndsi.dat %>% 
  arrange(recording.site, day, time, year) %>% 
  select(ndsi)
ordered.dat$ndsi <- ndsi.ordered$ndsi

# AEI
aei.dat <- extractMultipleSounds("../data/AEI.csv", index = 'AEI')

aei.ordered <- aei.dat %>% 
  arrange(recording.site, day, time, year) %>% 
  select(AEI)

ordered.dat$aei <- aei.ordered$AEI

# Treatment Indicator variable
ordered.dat$trt <- ifelse(ordered.dat$recording.site %in% c("LA09", "LA10", "LA11", "LA12"), 
                          1, 0)

# Read out full file for the ordered acoustic index data
write.csv(ordered.dat, "../data/orderedIndices.csv", row.names = FALSE)
