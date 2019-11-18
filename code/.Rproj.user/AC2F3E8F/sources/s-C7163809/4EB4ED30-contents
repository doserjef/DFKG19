# Clear working directory
rm(list =ls())
library(tuneR)
library(dplyr)
library(seewave)
library(soundecology)
# Contains a couple functions for manipulating output from soundecology functions
source("soundscape.R")

# This file computes the NDSI and the PSD for all recordings in a directory. 
# Other acoustic indices are computed separately

# Author: Jeffrey W. Doser

# Compute NDSI for all recordings -----------------------------------------

# Change directory as needed, as well as number of cores used
# Only need to run once, after run the first time can comment out
# multiple_sounds(directory = "../../soundsForJeff/ALL LA SITES-YEARS-0530-0700/", 
#                 resultfile = "../data/ndsi.csv",
#                 soundindex ="ndsi", no_cores = 4)

# Extract all necessary data  ---------------------------------------------

full.dat <- extractMultipleSounds("../data/ndsi.csv")
str(full.dat)
# Number of recordings
n <- dim(full.dat)[1]

# Compute PSD for all recordings ------------------------------------------

# Commented out to avoid long computations, only need to run once
# fileNames <- list.files(path = "../../soundsForJeff/ALL LA SITES-YEARS-0530-0700/", 
#                         pattern = "*.wav", full.names = TRUE)
# Compute PSDs for all recordings using psd function in the soundscape.R file
# psd.dat <- psd(fileNames)

# Write out the psd file so you only have to run once
# write.table(psd.dat, "../data/psd-raw", row.names = FALSE, col.names = FALSE, 
            # sep = "\t")

# Read in the psd file after initial computation
psd.dat <- read.table("../data/psd-raw")

# Combine psd and NDSI data
psd.dat <- cbind(psd.dat, full.dat)
  
# Arrange data by individual soundscape (i.e., by recordings site, day, time, year)
psd.ordered.dat <- psd.dat %>% 
  arrange(recording.site, day, time, year)

# Create a treatment indicator variable that takes value 1 for treatment sites and 
# 0 for control sites. 
psd.ordered.dat$trt <- ifelse(psd.ordered.dat$recording.site %in% 
                                c("LA09", "LA10", "LA11", "LA12"), 1, 0)

# Just grab the psd values 
psd.table <- psd.ordered.dat[, 1:10]
# Read out a table just containing psd values if desired
write.table(psd.table, "../data/psd-table.csv")
# Read out the ndsi and psd data in the properly ordered format
write.csv(psd.ordered.dat, "../data/full-ndsi-psd-ordered.csv")





