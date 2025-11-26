# This script runs a differential methylation analysis by genetic ancestry group across Horvath clock CpGs and creates a volcano plot of the results
# Author: Hao Xu & Rory Boyle rorytboyle@gmail.com
# Date: 24/11/2025

# Load required libraries
library(tidyverse)
library(readxl)
library(sesame)
library(Rtsne)
library(viridis)
library(parallel)
library(writexl)
library(circlize)
library(patchwork)
library(reshape2)
library(SummarizedExperiment)
library(limma)
library(ggrepel)
library(grid)
library(methylclockData) # assumes you have installed methylclockData package: BiocManager::install("methylclockData")


# read in Horvath weights from methylclockData package (verified that they are same coefficients as in Horvath 2013 Genome Biol paper here: HorvathClock_CpG_coefficients_sanity_check.R)
horvath_weights <- get_coefHorvath() %>%
  select(CpG = CpGmarker, Weights = CoefficientTraining) %>%
  # drop the intercept row
  filter(CpG != "(Intercept)")

# Load CpG sites of interest and filter data
cpg_file_path <- "/Users/rorytb/Library/CloudStorage/Box-Box/PMBB for Rory/DATA/20251124_cpg_Horvath_intersect_MSA.rds"
cpg_sites <- readRDS(cpg_file_path)

data <- readRDS("/Users/rorytb/Library/CloudStorage/Box-Box/PMBB for Rory/DATA/PMBBID_IDAT_covariate_class_05082025.rds")
data <- data %>% filter(Class %in% c("EUR", "AFR"))

# Load beta values and filter to match metadata
betasHM450_imputed <- readRDS('/Users/rorytb/Library/CloudStorage/Box-Box/PMBB for Rory/DATA/20241220_betasHM450_imputed.rds')
betasHM450_imputed <- betasHM450_imputed[, data$IDAT_file_name]  # Filter betas to match meta1 IDATs

# filter the betas to only include the sites used by the Horvath clock
cpg_sites <- unlist(cpg_sites)  # Convert list to vector
cpg_sites <- as.character(cpg_sites)  # Ensure it's a character vector
betas <- betasHM450_imputed[cpg_sites, , drop = FALSE]

# Differential methylation analysis using limma
group <- factor(data$Class, levels = c("EUR", "AFR"))
design <- model.matrix(~ group)
fit <- lmFit(betas, design)
fit <- eBayes(fit)

# Get all CpGs with stats
results <- topTable(fit, coef = 2, number = Inf, adjust.method = "fdr")

# Add delta beta (AFR - EUR)
results$delta_beta <- rowMeans(betas[, group == "AFR"], na.rm = TRUE) -
  rowMeans(betas[, group == "EUR"], na.rm = TRUE)

# Add Horvath weights to results
# Match weights by CpG names
cpg_names <- rownames(results)
results$horvath_weight <- horvath_weights$Weights[match(cpg_names, horvath_weights$CpG)]

# Check for NAs in weights
cat("NAs in horvath_weight:", sum(is.na(results$horvath_weight)), "\n")

# Create absolute weight for sizing
results$abs_weight <- abs(as.numeric(results$horvath_weight))

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

# Calculate plot limits
max_delta <- ceiling(max(abs(results$delta_beta * 100), na.rm = TRUE) / 5) * 5

# Create legend breaks that span the full range including 0
max_weight <- max(results$abs_weight, na.rm = TRUE)
legend_breaks <- c(0, 0.15, 0.3, 0.45, 0.6)

# Create volcano plot without annotations
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
  scale_color_manual(values = ancestry_colors) +
  guides(
    color = guide_legend(override.aes = list(size = 3, alpha = 1)),
    size = guide_legend(override.aes = list(color = "black", alpha = 1))
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  coord_cartesian(clip = "off") +
  annotate(
    "text",
    x = -Inf,
    y = -log10(0.05),
    label = "FDR = 0.05",
    hjust = 1.05,  # Right-align and push left of the axis
    vjust = 0.5,   # Center vertically on the line
    size = 14 / .pt,  # Match your fontsize
    color = "black"
  ) +
  scale_x_continuous(
    breaks = seq(-max_delta, max_delta, by = 5),
    limits = c(-max_delta-2, max_delta+2),  # Reduced from 8 to 2
    expand = expansion(mult = 0.02)  # Reduced from 0.05
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.05))  # Reduced top expansion
  ) +
  labs(
    x = expression(paste(Delta, "β(%)")),
    y = expression(-log[10](FDR)),
    color = "Methylation Level",
    size = "Clock Weight (absolute)"
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
    plot.margin = margin(5.5, 20, 5.5, 20, "pt")  # Reduced from 30 to 20
  )

volcano_by_ancestry

# Print summary statistics
n_total_sig <- sum(results_df$adj.P.Val < 0.05)
n_afr_sig <- sum(results_df$adj.P.Val < 0.05 & results_df$Ancestry == "Higher in AFR")
n_eur_sig <- sum(results_df$adj.P.Val < 0.05 & results_df$Ancestry == "Higher in EUR")

cat("Total significant CpGs:", n_total_sig, "\n")
cat("Significant CpGs higher in AFR:", n_afr_sig, "\n")
cat("Significant CpGs higher in EUR:", n_eur_sig, "\n")

# Save plot with narrower width
ggsave("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251124_Horvath_volcano_by_ancestry.png", 
       plot = volcano_by_ancestry, width = 8, height = 7, dpi = 300)  # Reduced from 11.4 to 8

# Save results
## add CpG as colname
results_df <- results_df %>%
  tibble::rownames_to_column(var = "CpG")

write.csv(results_df, "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251124_Horvath_DiffMethylAnalysis_results.csv", row.names = TRUE)