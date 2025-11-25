# This script summarizes the results of differential methylation analyses with and without adjusting for cell type proportions
# Author: Rory Boyle rorytboyle@gmail.com
# Date: 25/11/2025

# Load required libraries
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)

# Load the two results files
results_unadjusted <- read.csv("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251124_DunedinPACE_DiffMethylAnalysis_results.csv", row.names = 1)

results_adjusted <- read.csv("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251125_DunedinPACE_DiffMethylAnalysis_cell_type_prop_adjusted_results.csv", row.names = 1)

# Ensure CpG column is present and consistent
if (!"CpG" %in% colnames(results_unadjusted)) {
  results_unadjusted <- results_unadjusted %>%
    tibble::rownames_to_column(var = "CpG")
}

if (!"CpG" %in% colnames(results_adjusted)) {
  results_adjusted <- results_adjusted %>%
    tibble::rownames_to_column(var = "CpG")
}

# Merge the two datasets
comparison <- results_unadjusted %>%
  select(CpG, 
         logFC_unadj = logFC,
         AveExpr_unadj = AveExpr,
         t_unadj = t,
         P.Value_unadj = P.Value,
         adj.P.Val_unadj = adj.P.Val,
         B_unadj = B,
         delta_beta_unadj = delta_beta,
         dunedin_weight_unadj = dunedin_weight,
         abs_weight_unadj = abs_weight,
         Ancestry_unadj = Ancestry) %>%
  inner_join(
    results_adjusted %>%
      select(CpG,
             logFC_adj = logFC,
             AveExpr_adj = AveExpr,
             t_adj = t,
             P.Value_adj = P.Value,
             adj.P.Val_adj = adj.P.Val,
             B_adj = B,
             delta_beta_adj = delta_beta,
             dunedin_weight_adj = dunedin_weight,
             abs_weight_adj = abs_weight,
             Ancestry_adj = Ancestry),
    by = "CpG"
  )

# Calculate differences
comparison <- comparison %>%
  mutate(
    delta_logFC = logFC_adj - logFC_unadj,
    delta_t = t_adj - t_unadj,
    delta_delta_beta = (delta_beta_adj - delta_beta_unadj) * 100,  # Convert to %
    fold_change_pval = P.Value_adj / P.Value_unadj,
    significance_change = case_when(
      adj.P.Val_unadj < 0.05 & adj.P.Val_adj < 0.05 ~ "Sig in Both",
      adj.P.Val_unadj < 0.05 & adj.P.Val_adj >= 0.05 ~ "Lost Significance",
      adj.P.Val_unadj >= 0.05 & adj.P.Val_adj < 0.05 ~ "Gained Significance",
      TRUE ~ "Not Sig in Either"
    ),
    direction_change = case_when(
      sign(delta_beta_unadj) != sign(delta_beta_adj) ~ "Direction Changed",
      TRUE ~ "Direction Consistent"
    )
  )

# Summary statistics
print(table(comparison$significance_change))
print(table(comparison$direction_change))

# Plot Change in effect size (highlighting CpGs that were sig without adjusting for cell type prop and delta beta >5%)
comparison_labeled <- comparison %>%
  filter(
    adj.P.Val_unadj < 0.05 &  # FDR significant in unadjusted analysis
      abs(delta_beta_unadj * 100) > 5  # Delta beta > 5% in unadjusted
  )

print(comparison_labeled %>% select(CpG, delta_beta_unadj, adj.P.Val_unadj, delta_beta_adj, adj.P.Val_adj, delta_delta_beta))

# Create a combined status variable
comparison <- comparison %>%
  mutate(
    combined_status = case_when(
      adj.P.Val_unadj < 0.05 & adj.P.Val_adj < 0.05 & direction_change == "Direction Consistent" ~ "Sig in Both (Same Direction)",
      adj.P.Val_unadj < 0.05 & adj.P.Val_adj < 0.05 & direction_change == "Direction Changed" ~ "Sig in Both (Direction Changed)",
      adj.P.Val_unadj < 0.05 & adj.P.Val_adj >= 0.05 ~ "Lost Significance",
      adj.P.Val_unadj >= 0.05 & adj.P.Val_adj < 0.05 ~ "Gained Significance",
      TRUE ~ "Not Sig in Either"
    )
  )

# Update the labeled dataset
comparison_labeled <- comparison %>%
  filter(
    adj.P.Val_unadj < 0.05 &
      abs(delta_beta_unadj * 100) > 5
  )

delta_beta_change_plot <- ggplot(comparison, aes(x = delta_beta_unadj * 100, 
                             y = delta_delta_beta,
                             color = combined_status)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.6, size = 2) +
  geom_text_repel(
    data = comparison_labeled,
    aes(label = CpG),
    size = 6,
    max.overlaps = 30,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = 'grey50',
    min.segment.length = 0
  ) +
  scale_color_manual(
    values = c(
      "Sig in Both (Same Direction)" = "#0072B2",
      "Sig in Both (Direction Changed)" = "#D55E00",
      "Lost Significance" = "#CC79A7",
      "Gained Significance" = "#F0E442",
      "Not Sig in Either" = "grey70"
    )
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 4, alpha = 1))
  ) +
  labs(
    x = expression(paste(Delta, "β(%) - Unadjusted")),
    y = expression(paste("Change in ", Delta, "β(%) after cell type adjustment")),
    color = "Status after adjustment\nfor cell type proportions",
    title = "Δβ Change of DunedinPACE CpGs by genetic ancestry after cell type adjustment",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.subtitle = element_text(size = 10, color = "grey40")
  )

delta_beta_change_plot

ggsave(
  filename = "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251125_DunedinPACE_DiffMethylAnalysis_cell_type_adjustment_delta_beta_change_plot.png",
  plot = delta_beta_change_plot,
  width = 8,
  height = 6,
  dpi = 300
)
