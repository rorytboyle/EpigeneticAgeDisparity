# This script runs a sanity check of Horvath (353) clock weights to ensure that the same model coefficients are received 
# from the methylclock github, methylclockdata package and the supplementary data provided in the Horvath 2013 Genome Biol paper: https://doi.org/10.1186/gb-2013-14-10-r115

# Author: Rory Boyle rorytboyle@gmail.com
# Date: 24/11/2025

library(tidyverse)
library(methylclockData) # assumes you have installed methylclockData package: BiocManager::install("methylclockData")

# Download Horvath clock coefficients from methylclock github
download.file("https://github.com/isglobal-brge/methylclock/raw/main/data/coefHorvath.rda",
              destfile = "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/data/coefHorvath.rda",
              mode = "wb")

load("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/data/coefHorvath.rda")
github_coeffs <- coefHorvath

# Read in coeffs from Horvath 2013 manuscript
manu_coeffs <- read.csv("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/data/13059_2013_3156_MOESM3_ESM.csv",
                        skip = 2, header = TRUE)

# Read in coeffs from methylclockData bioconductor package
methylclockData_coeffs <- get_coefHorvath()

# verify that coefficient weights are the same in each
all(github_coeffs$CoefficientTraining == manu_coeffs$CoefficientTraining)
all(github_coeffs$CoefficientTraining == methylclockData_coeffs$CoefficientTraining)

#verify the CpGs are the asme in each
all(github_coeffs$CpGmarker == manu_coeffs$CpGmarker)
all(github_coeffs$CpGmarker == methylclockData_coeffs$CpG)
