# This script runs a cell type deconvolution analysis on PMBB DNA methylation data using the EpiDISH package
# Author: Hao Xu
# Date: 25/11/2025

library(EpiDISH) # assumes you have installed methylclockData package: BiocManager::install("EpiDISH")
library(tidyverse)

# Load reference data for blood cell types
data(centDHSbloodDMC.m)

# Load metadata for DNAm data
data <- readRDS("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/data/PMBBID_IDAT_covariate_class_05082025.rds")

# Load beta values and filter to match metadata
betasHM450_imputed <- readRDS("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/data/20241220_betasHM450_imputed.rds")
dim(betasHM450_imputed)  # Check dimensions
# [1] 486427    848
betasHM450_imputed <- betasHM450_imputed[, data$IDAT_file_name]  # Filter betas to match meta1 IDATs
dim(betasHM450_imputed)  # Check dimensions
# [1] 486427    811

# Run EpiDISH to estimate cell type proportions
BloodFrac.m <- epidish(beta.m = betasHM450_imputed, ref.m = centDHSbloodDMC.m, method = "RPC")$estF

# Save the cell type proportions
saveRDS(BloodFrac.m, file="/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/data/20250916_BloodFrac_m.rds")

# save the boxplot
BloodFrac_df <- as.data.frame(BloodFrac.m)
BloodFrac_df$SampleID <- rownames(BloodFrac_df)
BloodFrac_long <- reshape2::melt(BloodFrac_df, id.vars = "SampleID", variable.name = "CellType", value.name = "Proportion")
p <- ggplot(BloodFrac_long, aes(x = CellType, y = Proportion)) +
    geom_boxplot() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Cell Type Proportions (811 samples)", x = "Cell Type", y = "Proportion")

ggsave("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251125_BloodFrac_boxplot.pdf", plot = p, width = 8, height = 6)
