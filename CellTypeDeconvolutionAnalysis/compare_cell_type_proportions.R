library(ggplot2)
library(reshape2)
library(ggpubr)
library(sjPlot)
cbPalette <- c(
  "#E69F00",  # orange
  "#009E73",  # bluish green
  "#7F3C3C",  # dark red (custom)
  "#56B4E9",  # sky blue
  "#CC79A7",  # reddish purple
  "#F0E442"   # yellow
)

# Prep data ####
# Read in epidish output
epidish_df <- readRDS('/Users/rorytb/Library/CloudStorage/Box-Box/PMBB for Rory/cell_type_deconvolution/20250915_BloodFrac_m.rds') %>%
  as.data.frame() %>%
  tibble::rownames_to_column("IDAT_file_name")

# Read in IDAT file names and covariate info
idat_files <- readRDS('/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/PMBBID_IDAT_covariate_class_05082025.rds') %>%
  select(PMBB_ID, IDAT_file_name, Class)

# Merge idat_files with long_data
epidish_df <- merge(epidish_df, idat_files, by = "IDAT_file_name")

# Statistical analysis ####
# Specify cell type columns
cell_type_cols <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino")

# Specify your two groups here
group1_name <- "AFR"  # Replace with your actual group name
group2_name <- "EUR"  # Replace with your actual group name

# Create results data frame
results <- data.frame(
  cell_type = cell_type_cols,
  p_value = NA,
  mean_group1 = NA,
  mean_group2 = NA,
  test_method = NA,
  effect_size = NA
)

# Statistical testing for each cell type
for(i in 1:length(cell_type_cols)) {
  cell_type <- cell_type_cols[i]
  
  # Extract proportions for each group
  group1_props <- epidish_df[epidish_df$Class == group1_name, cell_type]
  group2_props <- epidish_df[epidish_df$Class == group2_name, cell_type]
  
  # Remove any NA values
  group1_props <- group1_props[!is.na(group1_props)]
  group2_props <- group2_props[!is.na(group2_props)]
  
  # Check normality and choose appropriate test
  if(length(group1_props) > 3 && length(group2_props) > 3) {
    norm_test1 <- shapiro.test(group1_props)$p.value > 0.05
    norm_test2 <- shapiro.test(group2_props)$p.value > 0.05
    
    if(norm_test1 && norm_test2) {
      # Use t-test if normally distributed
      test_result <- t.test(group1_props, group2_props)
      results$test_method[i] <- "t-test"
    } else {
      # Use Wilcoxon test if not normally distributed
      test_result <- wilcox.test(group1_props, group2_props)
      results$test_method[i] <- "Wilcoxon"
    }
  } else {
    # Use Wilcoxon for small samples
    test_result <- wilcox.test(group1_props, group2_props)
    results$test_method[i] <- "Wilcoxon"
  }
  
  # Store results
  results$p_value[i] <- test_result$p.value
  results$mean_group1[i] <- mean(group1_props)
  results$mean_group2[i] <- mean(group2_props)
  
  # Calculate effect size (Cohen's d)
  pooled_sd <- sqrt(((length(group1_props)-1)*var(group1_props) + 
                       (length(group2_props)-1)*var(group2_props)) / 
                      (length(group1_props) + length(group2_props) - 2))
  cohens_d <- (mean(group1_props) - mean(group2_props)) / pooled_sd
  results$effect_size[i] <- cohens_d
}

# Apply multiple testing correction
results$p_adjusted <- p.adjust(results$p_value, method = "fdr")

# Add group names to results for clarity
results$group1 <- "AFR"
results$group2 <- "EUR"

print(results)

# Visualisation ####


# Filter data to only include your two specified groups
epidish_filtered <- epidish_df[epidish_df$Class %in% c(group1_name, group2_name), ]

# Create long format for plotting
plot_data <- melt(epidish_filtered[, c("Class", cell_type_cols)], 
                  id.vars = "Class",
                  variable.name = "Cell_Type", 
                  value.name = "Proportion")

# Create violin plots as alternative
cell_type_ancestry_plot <- ggboxplot(
  plot_data, 
  x = "Class", 
  y = "Proportion",
  color = "Class",
  palette = cbPalette,
  add = "jitter",
  ggtheme = theme_sjplot2()) +
  facet_wrap(~Cell_Type, scales = "free_y", ncol = 3) +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  ) +
  labs(
    y = "Proportion"
  ) +
  stat_compare_means(
    label = "p.format",  
    label.x.npc = "center",
    vjust = 1,    # Move labels further up
    size = 4,      # Larger text
    label.y.npc = 0.95  # Position near top of each panel
  )
    
cell_type_ancestry_plot
