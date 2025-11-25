# This script runs a differential methylation analysis by genetic ancestry group adjusted by cell type proportion across DunedinPACE CpGs and creates a volcano plot of the results
# Author: Hao Xu & Rory Boyle rorytboyle@gmail.com
# Date: 24/11/2025

# Load required libraries
library(dplyr)
library(ggplot2)
library(limma)
library(ggrepel)
library(grid)
library(DunedinPACE)

cell_prop <- readRDS("/home/xuh6/20240920_Hao/projects/PMBB_data_analysis/20250916_BloodFrac_m.rds")

dunedinPACE_weights <- mPACE_Models$model_weights

# Load CpG sites of interest and filter data
cpg_file_path <- "/home/xuh6/20240920_Hao/projects/PMBB_data_analysis/20250529_cpg_DunedinPACE_intersect_MSA.rds"
cpg_sites <- readRDS(cpg_file_path)

data <- readRDS("/home/xuh6/20240920_Hao/projects/PMBB_data_analysis/PMBBID_IDAT_covariate_class_05082025.rds")
data <- data %>% filter(Class %in% c("EUR", "AFR"))

# Load beta values and filter to match metadata
betasHM450_imputed <- readRDS("/home/xuh6/20240920_Hao/projects/PMBB_data_analysis/20241220_betasHM450_imputed.rds")
betasHM450_imputed <- betasHM450_imputed[, data$IDAT_file_name]

# Make sure cell_prop match the order in beta matrix and metadata
cell_prop <- cell_prop[data$IDAT_file_name, , drop = FALSE]
# Check alignment (should return TRUE)
all(rownames(cell_prop) == data$IDAT_file_name)

# Filter betas to only include the sites used by the DunedinPACE clock
cpg_sites <- unlist(cpg_sites)
cpg_sites <- as.character(cpg_sites)
betas <- betasHM450_imputed[cpg_sites, , drop = FALSE]

# Differential methylation analysis using limma
group <- factor(data$Class, levels = c("EUR", "AFR"))
design <- model.matrix(~ group + as.matrix(cell_prop))
# remove the col as.matrix(cell_prop)Eosino
design <- design[, !grepl("Eosino", colnames(design))]
fit <- lmFit(betas, design)
fit <- eBayes(fit)

# Get all CpGs with stats
results <- topTable(fit, coef = 2, number = Inf, adjust.method = "fdr")

# Add delta beta (AFR - EUR)
results$delta_beta <- rowMeans(betas[, group == "AFR"], na.rm = TRUE) -
  rowMeans(betas[, group == "EUR"], na.rm = TRUE)

# Add DunedinPACE weights
cpg_names <- rownames(results)
results$dunedin_weight <- dunedinPACE_weights$DunedinPACE[cpg_names]

# Check for NAs in weights
cat("NAs in dunedin_weight:", sum(is.na(results$dunedin_weight)), "\n")

# Create absolute weight for sizing
results$abs_weight <- abs(as.numeric(results$dunedin_weight))

# Classify by ancestry
results <- results %>%
  mutate(Ancestry = case_when(
    adj.P.Val >= 0.05 ~ "Not Significant",
    delta_beta > 0 ~ "Higher in AFR",
    delta_beta < 0 ~ "Higher in EUR",
    TRUE ~ "Not Significant"
  )) %>%
  mutate(Ancestry = factor(Ancestry, levels = c("Higher in AFR", "Higher in EUR", "Not Significant")))

# Define color palette
ancestry_colors <- c(
  "Higher in AFR" = "#E69F00",  # orange
  "Higher in EUR" = "#009E73",  # bluish green
  "Not Significant" = "grey"
)

# Convert to data.frame
results_df <- as.data.frame(results)

# Filter for significant points to label
results_sig <- results_df %>%
  filter(adj.P.Val < 0.05, abs(delta_beta * 100) > 5) 

# Calculate plot limits
max_delta <- ceiling(max(abs(results$delta_beta * 100), na.rm = TRUE) / 5) * 5

# Create legend breaks that span the full range including 0
max_weight <- max(results$abs_weight, na.rm = TRUE)
legend_breaks <- c(0, 0.15, 0.3, 0.45, 0.6)

# Create volcano plot
volcano_by_ancestry <- ggplot(results, aes(x = delta_beta * 100, 
                                           y = -log10(adj.P.Val), 
                                           color = Ancestry, 
                                           size = abs_weight)) +
  geom_point(alpha = 0.6, shape = 16) +
  scale_size_continuous(
    range = c(2, 6),  # min and max point sizes
    breaks = legend_breaks,
    labels = function(x) sprintf("%.2f", x),
    limits = c(0, NA)  # Ensure scale starts at 0
  ) +
  geom_text_repel(
    data = results_sig,
    aes(label = rownames(results_sig)),
    size = 5,
    max.overlaps = 20,
    box.padding = 0.4,
    point.padding = 0.5,
    segment.color = 'black',
    nudge_x = case_when(
      rownames(results_sig) == "cg00151250" ~ 4,
      results_sig$delta_beta > 0 ~ 3,
      TRUE ~ -4
    ),
    nudge_y = ifelse(rownames(results_sig) == "cg00151250", 7, 0.5)
  ) +
  scale_color_manual(values = ancestry_colors) +
  guides(
    color = guide_legend(override.aes = list(size = 3, alpha = 1)),
    size = guide_legend(override.aes = list(color = "black", alpha = 1))
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  annotation_custom(
    grob = grid::textGrob(
      label = "FDR = 0.05", 
      x = unit(-0.19, "npc"),
      y = unit(0.03, "npc"),
      just = "left",
      gp = gpar(fontsize = 14)
    ),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  coord_cartesian(clip = "off") +
  geom_vline(xintercept = c(-5, 5), linetype = "dashed", color = "grey") +
  scale_x_continuous(
    breaks = seq(-max_delta, max_delta, by = 5),
    limits = c(-max_delta-5, max_delta+5),
    expand = expansion(mult = 0.02)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = expression(paste(Delta, "β(%)")),
    y = expression(-log[10](FDR)),
    color = "Methylation Level",
    size = "DunedinPACE\nWeight (absolute)"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = unit(0.2, "cm"),
    plot.margin = margin(5.5, 20, 5.5, 20, "pt")
  )

volcano_by_ancestry
# Print summary statistics
n_total_sig <- sum(results_df$adj.P.Val < 0.05)
n_afr_sig <- sum(results_df$adj.P.Val < 0.05 & results_df$Ancestry == "Higher in AFR")
n_eur_sig <- sum(results_df$adj.P.Val < 0.05 & results_df$Ancestry == "Higher in EUR")

cat("Total significant CpGs:", n_total_sig, "\n")
cat("Significant CpGs higher in AFR:", n_afr_sig, "\n")
cat("Significant CpGs higher in EUR:", n_eur_sig, "\n")

# Save plot
ggsave("20251124_171_volcano_by_ancestry_adjusted_cell_type.png", 
       plot = volcano_by_ancestry, width = 10, height = 6, dpi = 300)

# Save results
## add CpG as colname
results_df <- results_df %>%
  tibble::rownames_to_column(var = "CpG")

write.csv(results_df, "20251124_171_DiffMethylAnalysis_results_adjusted_cell_type.csv", row.names = TRUE)
