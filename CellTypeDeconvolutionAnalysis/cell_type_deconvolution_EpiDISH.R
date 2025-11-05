library(EpiDISH)
library(tidyr)
library(dplyr)
data(centDHSbloodDMC.m)

data <- readRDS("/home/xuh6/20240920_Hao/projects/PMBB_data_analysis/PMBBID_IDAT_covariate_class_05082025.rds")

# Load beta values and filter to match metadata
betasHM450_imputed <- readRDS("20241220_betasHM450_imputed.rds")
dim(betasHM450_imputed)  # Check dimensions
# [1] 486427    848
betasHM450_imputed <- betasHM450_imputed[, data$IDAT_file_name]  # Filter betas to match meta1 IDATs
dim(betasHM450_imputed)  # Check dimensions
# [1] 486427    811
BloodFrac.m <- epidish(beta.m = betasHM450_imputed, ref.m = centDHSbloodDMC.m, method = "RPC")$estF

# Save the cell type proportions
saveRDS(BloodFrac.m, file="20250916_BloodFrac_m.rds")

# save the boxplot
BloodFrac_df <- as.data.frame(BloodFrac.m)
BloodFrac_df$SampleID <- rownames(BloodFrac_df)
BloodFrac_long <- reshape2::melt(BloodFrac_df, id.vars = "SampleID", variable.name = "CellType", value.name = "Proportion")
library(ggplot2)
p <- ggplot(BloodFrac_long, aes(x = CellType, y = Proportion)) +
    geom_boxplot() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Cell Type Proportions (811 samples)", x = "Cell Type", y = "Proportion")
ggsave("20250916_BloodFrac_boxplot.pdf", plot = p, width = 8, height = 6)
