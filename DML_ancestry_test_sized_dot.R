# Load required libraries
library(readxl)
library(sesame)
library(dplyr)
library(Rtsne)
library(ggplot2)
library(viridis)
library(parallel)
library(writexl)
library(circlize)
library(patchwork)
library(reshape2)
library(SummarizedExperiment)
library(tidyr)
library(limma)
library(ggrepel)
library(grid)
library(DunedinPACE)

dunedinPACE_weights <- mPACE_Models$model_weights


# Load CpG sites of interest and filter data
cpg_file_path <- "/home/xuh6/20240920_Hao/projects/PMBB_data_analysis/20250529_cpg_DunedinPACE_intersect_MSA.rds"
cpg_sites <- readRDS(cpg_file_path)

data <- readRDS("/home/xuh6/20240920_Hao/projects/PMBB_data_analysis/PMBBID_IDAT_covariate_class_05082025.rds")
data <- data %>% filter(Class %in% c("EUR", "AFR"))

# Load beta values and filter to match metadata
betasHM450_imputed <- readRDS("20241220_betasHM450_imputed.rds")
betasHM450_imputed <- betasHM450_imputed[, data$IDAT_file_name]  # Filter betas to match meta1 IDATs


# Load the CpG sites and betas
cpg_sites <- readRDS(cpg_file_path)
# filter the betas to only include the sites used by the DunedinPACE clock
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

# Convert delta beta to percent
results$delta_beta_pct <- results$delta_beta * 100

# Add DunedinPACE weights to results
# Match weights by CpG names
cpg_names <- rownames(results)
results$dunedin_weight <- dunedinPACE_weights$DunedinPACE[cpg_names]

# Handle missing weights (set to 0 or minimum weight)
# check if there are any NAs in dunedin_weight
print("check NAs in dunedin_weight:")
sum(is.na(results$dunedin_weight))
# results$dunedin_weight[is.na(results$dunedin_weight)] <- 0

# Create absolute weight for sizing (since weights can be negative)
results$abs_weight <- abs(as.numeric(results$dunedin_weight))

# Scale the absolute weights for better visualization
min_size <- 2        # Minimum point size
max_size <- 6           # Maximum point size

# Scale weights to size range
if(max(results$abs_weight, na.rm = TRUE) > 0) {
  results$point_size <- min_size + (results$abs_weight / max(results$abs_weight, na.rm = TRUE)) * (max_size - min_size)
} else {
  results$point_size <- min_size
}

# Define thresholds
delta_thresh <- 0.05     # 20% methylation difference
fdr_thresh <- 0.05

results <- results %>%
  mutate(Ancestry = case_when(
    adj.P.Val >= 0.05 ~ "Not Significant",
    delta_beta > 0 ~ "Higher in AFR",
    delta_beta < -0 ~ "Higher in EUR",
    TRUE ~ "Not Significant"
  )) %>%
  mutate(Ancestry = factor(Ancestry, levels = c("Higher in AFR", "Higher in EUR", "Not Significant")))


# Define new color palette
ancestry_colors <- c(
  "Higher in AFR" = "orange",
  "Higher in EUR" = "green",
  "Not Significant" = "grey"
)
# Convert to proper data.frame
results_df <- as.data.frame(results)

export_df <- results_df %>%
  mutate(CpG = rownames(results_df)) %>%
  select(CpG, P.Value, adj.P.Val, delta_beta, Ancestry)
# write.csv(export_df, "20250819_171_stats.csv", row.names = FALSE)

# Filter for significant points and sort
results_sig <- results_df %>%
  filter(adj.P.Val < 0.05, abs(delta_beta * 100) > 5) 

max_delta <- ceiling(max(abs(results$delta_beta * 100), na.rm = TRUE) / 5) * 5
# Plot
volcano_by_ancestry <- ggplot(results, aes(x = delta_beta * 100, y = -log10(adj.P.Val), color = Ancestry,  size = point_size)) +
  geom_point(alpha = 0.6, shape = 16) +
  scale_size_identity() +
    geom_text_repel(
      data = results_sig,
      aes(label = rownames(results_sig)),
      size = 5,
      max.overlaps = 20,
      box.padding = 0.4,
      point.padding = 0.5,
      segment.color = 'black',
      nudge_x = case_when(
        rownames(results_sig) == "cg00151250" ~ 4,               # move right
        # rownames(results_sig) == "cg25243766" ~ 8,               # move right
        results_sig$delta_beta > 0 ~ 3,
        TRUE ~ -4
      ),
      nudge_y = ifelse(rownames(results_sig) == "cg00151250", 7, 0.5)  # push cg00151250 higher
  ) +
  scale_color_manual(values = ancestry_colors) +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 1.5))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  annotation_custom(
    grob = grid::textGrob(
      label = "FDR = 0.05", 
      x = unit(-0.19, "npc"),  # absolute left of plot panel
      y = unit(0.03, "npc"),  # around y = -log10(0.05)
      just = "left",
      gp = gpar(fontsize = 14)
    ),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  coord_cartesian(clip = "off") +
  geom_vline(xintercept = delta_thresh*100, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = -delta_thresh*100, linetype = "dashed", color = "grey") +
  scale_x_continuous(
    breaks = seq(-max_delta, max_delta, by = 5),
    limits = c(-max_delta-5, max_delta+5),
    expand = expansion(mult = 0.02),
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = expression(paste(Delta, "β(%)")),
    y = expression(-log[10](FDR)),
    color = "Methylation Level"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 14),       # Size of legend labels
    legend.title = element_text(size = 16),       # Size of legend title
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.ticks.length = unit(0.2, "cm"),
    plot.margin = margin(5.5, 5.5, 5.5, 40, "pt")
  )

# Display
print(volcano_by_ancestry)

# Optional: Save plot
ggsave("20250915_171_volcano_by_ancestry.png", plot = volcano_by_ancestry, width = 10, height = 6, dpi = 300)

n_total_sig <- results_df %>%
  filter(adj.P.Val < 0.05) %>%
  nrow()
n_afr_sig <- results_df %>%
  filter(adj.P.Val < 0.05, Ancestry == "Higher in AFR") %>%
  nrow()
n_eur_sig <- results_df %>%
  filter(adj.P.Val < 0.05, Ancestry == "Higher in EUR") %>%
  nrow()
cat("Total significant CpGs:", n_total_sig, "\n")
cat("Significant CpGs higher in AFR:", n_afr_sig, "\n")
cat("Significant CpGs higher in EUR:", n_eur_sig, "\n")








# Differential methylation analysis using limma
group <- factor(data$Class, levels = c("EUR", "AFR"))
design <- model.matrix(~ group)
fit <- lmFit(betas, design)
fit <- eBayes(fit)

results <- topTable(fit, coef = 2, number = Inf, adjust.method = "none")

# Add delta beta (AFR - EUR)
results$delta_beta <- rowMeans(betas[, group == "AFR"], na.rm = TRUE) -
                      rowMeans(betas[, group == "EUR"], na.rm = TRUE)

# Convert delta beta to percent
results$delta_beta_pct <- results$delta_beta * 100

# Add DunedinPACE weights to results
# Match weights by CpG names
cpg_names <- rownames(results)
results$dunedin_weight <- dunedinPACE_weights$DunedinPACE[cpg_names]

# Handle missing weights (set to 0 or minimum weight)
# check if there are any NAs in dunedin_weight
print("check NAs in dunedin_weight:")
sum(is.na(results$dunedin_weight))

# Create absolute weight for sizing (since weights can be negative)
results$abs_weight <- abs(as.numeric(results$dunedin_weight))

# Scale the absolute weights for better visualization
min_size <- 2        # Minimum point size
max_size <- 6           # Maximum point size

# Scale weights to size range
if(max(results$abs_weight, na.rm = TRUE) > 0) {
  results$point_size <- min_size + (results$abs_weight / max(results$abs_weight, na.rm = TRUE)) * (max_size - min_size)
} else {
  results$point_size <- min_size
}

results <- results %>%
  mutate(Ancestry = case_when(
    adj.P.Val >= 0.05 ~ "Not Significant",
    delta_beta > 0 ~ "Higher in AFR",
    delta_beta < -0 ~ "Higher in EUR",
    TRUE ~ "Not Significant"
  )) %>%
  mutate(Ancestry = factor(Ancestry, levels = c("Higher in AFR", "Higher in EUR", "Not Significant")))


# Define new color palette
ancestry_colors <- c(
  "Higher in AFR" = "orange",
  "Higher in EUR" = "green",
  "Not Significant" = "grey"
)

# Convert to proper data.frame
results_df <- as.data.frame(results)

# Filter for significant points and sort
results_sig <- results_df %>%
  filter(adj.P.Val < 0.05, abs(delta_beta * 100) > 5) 

max_delta <- ceiling(max(abs(results$delta_beta * 100), na.rm = TRUE) / 5) * 5
# Plot
volcano_by_ancestry <- ggplot(results, aes(x = delta_beta * 100, y = -log10(adj.P.Val), color = Ancestry,  size = point_size)) +
  geom_point(alpha = 0.6, shape = 16) +
  scale_size_identity() +
    geom_text_repel(
      data = results_sig,
      aes(label = rownames(results_sig)),
      size = 5,
      max.overlaps = 20,
      box.padding = 0.4,
      point.padding = 0.5,
      segment.color = 'black',
      # nudge_x = ifelse(results_sig$delta_beta > 0, 2, -2),  # smart nudging left/right
      nudge_x = case_when(
        rownames(results_sig) == "cg00151250" ~ 4,               # move right
        # rownames(results_sig) == "cg25243766" ~ 8,               # move right
        results_sig$delta_beta > 0 ~ 3,
        TRUE ~ -4
      ),
      nudge_y = ifelse(rownames(results_sig) == "cg00151250", 7, 0.5)  # push cg00151250 higher
  ) +
  scale_color_manual(values = ancestry_colors) +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 1.5))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  annotation_custom(
    grob = grid::textGrob(
      label = "p = 0.05", 
      x = unit(-0.145, "npc"),  # absolute left of plot panel
      y = unit(0.03, "npc"),  # around y = -log10(0.05)
      just = "left",
      gp = gpar(fontsize = 14)
    ),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  coord_cartesian(clip = "off") +
  geom_vline(xintercept = delta_thresh*100, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = -delta_thresh*100, linetype = "dashed", color = "grey") +
  scale_x_continuous(
    breaks = seq(-max_delta, max_delta, by = 5),
    limits = c(-max_delta-5, max_delta+5),
    expand = expansion(mult = 0.02),
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    # title = "Volcano Plot of CpG Methylation by Ancestry",
    x = expression(paste(Delta, "β(%)")),
    y = expression(-log[10](p)),
    color = "Methylation Level"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 14),       # Size of legend labels
    legend.title = element_text(size = 16),       # Size of legend title
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.ticks.length = unit(0.2, "cm"),
    plot.margin = margin(5.5, 5.5, 5.5, 40, "pt")
  )

# Display
print(volcano_by_ancestry)

# Optional: Save plot
ggsave("20250915_171_volcano_by_ancestry_unadjusted.png", plot = volcano_by_ancestry, width = 10, height = 6, dpi = 300)

n_total_sig <- results_df %>%
  filter(adj.P.Val < 0.05) %>%
  nrow()
n_afr_sig <- results_df %>%
  filter(adj.P.Val < 0.05, Ancestry == "Higher in AFR") %>%
  nrow()
n_eur_sig <- results_df %>%
  filter(adj.P.Val < 0.05, Ancestry == "Higher in EUR") %>%
  nrow()
cat("Total significant CpGs:", n_total_sig, "\n")
cat("Significant CpGs higher in AFR:", n_afr_sig, "\n")
cat("Significant CpGs higher in EUR:", n_eur_sig, "\n")



