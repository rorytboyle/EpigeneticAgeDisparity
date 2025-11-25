# This script summarizes the results of differential methylation analyses with and without adjusting for cell type proportions
# Author: Rory Boyle rorytboyle@gmail.com
# Date: 25/11/2025

# Load required libraries
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)

# Function to compare unadjusted DML results vs adjusted for cell type proportions ####
compare_adjusted_unadjusted <- function(unadj_file, adj_file, clock_name, weight_col) {
  
  cat(paste("ANALYZING", toupper(clock_name), "CLOCK\n"))
  cat("\n")
  
  # Load the two results files
  results_unadjusted <- read.csv(unadj_file, row.names = 1)
  results_adjusted <- read.csv(adj_file, row.names = 1)
  
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
           weight_unadj = !!sym(weight_col),
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
               weight_adj = !!sym(weight_col),
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
  cat("Significance changes:\n")
  print(table(comparison$significance_change))
  cat("\n")
  
  cat("Direction changes:\n")
  print(table(comparison$direction_change))
  cat("\n")
  
  # Create a combined status variable (for plotting)
  comparison <- comparison %>%
    mutate(
      combined_status = case_when(
        adj.P.Val_unadj < 0.05 & adj.P.Val_adj < 0.05 & direction_change == "Direction Consistent" ~ "Sig in Both (Same Direction)",
        adj.P.Val_unadj < 0.05 & adj.P.Val_adj < 0.05 & direction_change == "Direction Changed" ~ "Sig in Both (Direction Changed)",
        adj.P.Val_unadj < 0.05 & adj.P.Val_adj >= 0.05 ~ "Lost Significance",
        adj.P.Val_unadj >= 0.05 & adj.P.Val_adj < 0.05 ~ "Gained Significance",
        TRUE ~ "Not Sig in Either"
      ),
      # Create a cleaner shape variable for direction - removed (unadj) suffix
      direction_unadj = case_when(
        Ancestry_unadj == "Higher in AFR" ~ "Higher in AFR",
        Ancestry_unadj == "Higher in EUR" ~ "Higher in EUR",
        TRUE ~ "Not Significant"
      )
    )
  
  # Filter for CpGs that were sig and beyond 5% delta beta in unadjusted differential methylation analyses
  comparison_labeled <- comparison %>%
    filter(
      adj.P.Val_unadj < 0.05 &
        abs(delta_beta_unadj * 100) > 5
    )
  
  cat("Number of CpGs to be labeled:", nrow(comparison_labeled), "\n")
  cat("Labeled CpGs:\n")
  print(comparison_labeled %>% select(CpG, delta_beta_unadj, adj.P.Val_unadj, delta_beta_adj, adj.P.Val_adj, delta_delta_beta))
  cat("\n")
  
  # Create plot with both color and shape
  delta_beta_change_plot <- ggplot(comparison, aes(x = delta_beta_unadj * 100, 
                                                   y = delta_delta_beta,
                                                   color = combined_status,
                                                   shape = direction_unadj)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.6, size = 2.5) +
    geom_text_repel(
      data = comparison_labeled,
      aes(label = CpG),
      size = 4,
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
    scale_shape_manual(
      values = c(
        "Higher in AFR" = 16,  # Filled circle
        "Higher in EUR" = 17,  # Filled triangle
        "Not Significant" = 15  # Filled square
      )
    ) +
    guides(
      shape = guide_legend(override.aes = list(size = 4, alpha = 1), order = 1),
      color = guide_legend(override.aes = list(size = 4, alpha = 1), order = 2)
    ) +
    labs(
      x = expression(paste(Delta, "β(%) - Unadjusted")),
      y = expression(paste("Change in ", Delta, "β(%) after cell type adjustment")),
      shape = "Direction without adjusting\nfor cell type proportions",
      color = "Status after adjustment\nfor cell type proportions",
      title = paste0("Δβ Change of ", clock_name, " CpGs by genetic ancestry after cell type adjustment")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      legend.box = "vertical"
    )
  
  return(list(
    comparison = comparison,
    plot = delta_beta_change_plot,
    clock_name = clock_name
  ))
}

# Run comparisons for both clocks ####
# DunedinPACE
dunedinpace_results <- compare_adjusted_unadjusted(
  unadj_file = "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251124_DunedinPACE_DiffMethylAnalysis_results.csv",
  adj_file = "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251125_DunedinPACE_DiffMethylAnalysis_cell_type_prop_adjusted_results.csv",
  clock_name = "DunedinPACE",
  weight_col = "dunedin_weight"
)

# Print DunedinPACE plot
print(dunedinpace_results$plot)

# Horvath
horvath_results <- compare_adjusted_unadjusted(
  unadj_file = "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251124_Horvath_DiffMethylAnalysis_results.csv",
  adj_file = "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251125_Horvath_DiffMethylAnalysis_cell_type_prop_adjusted_results.csv",
  clock_name = "Horvath",
  weight_col = "horvath_weight"
)

# Print Horvath plot
print(horvath_results$plot)

# Save plots ####
# DunedinPACE
ggsave(
  filename = "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251125_DunedinPACE_DiffMethylAnalysis_cell_type_adjustment_delta_beta_change_plot.png",
  plot = dunedinpace_results$plot,
  width = 11,
  height = 7,
  dpi = 300
)

# Horvath
ggsave(
  filename = "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251125_Horvath_DiffMethylAnalysis_cell_type_adjustment_delta_beta_change_plot.png",
  plot = horvath_results$plot,
  width = 11,
  height = 7,
  dpi = 300
)

# Combined side-by-side plot with single legend
combined_plot <- (dunedinpace_results$plot + 
                    theme(legend.position = "none") + 
                    labs(title = "DunedinPACE")) + 
  (horvath_results$plot + 
     labs(title = "Horvath")) +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Δβ Change of Clock CpGs by Genetic Ancestry After Cell Type Adjustment",
    theme = theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))
  )

combined_plot

ggsave(
  filename = "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251125_Combined_DiffMethylAnalysis_cell_type_adjustment_delta_beta_change_plot.png",
  plot = combined_plot,
  width = 18,
  height = 7,
  dpi = 300
)

# Save comparison data ####

# DunedinPACE
write.csv(
  dunedinpace_results$comparison,
  "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251125_DunedinPACE_adjusted_vs_unadjusted_comparison.csv",
  row.names = FALSE
)

# Horvath
write.csv(
  horvath_results$comparison,
  "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251125_Horvath_adjusted_vs_unadjusted_comparison.csv",
  row.names = FALSE
)