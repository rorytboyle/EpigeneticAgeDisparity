# This script creates a volcano plot of a differential methylation analysis that
# annotates hypermethylated CpGs that were enriched for genes associated with neurodegenerative traits.
# Author: Rory Boyle rorytboyle@gmail.com
# Date: 09/12/2025

# Load required libraries
library(tidyverse)
library(ggrepel)
library(sesame)
library(sesameData)
library(knowYourCG)

# Load the differential methylation results
results <- read.csv("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251125_Horvath_DiffMethylAnalysis_cell_type_prop_adjusted_results.csv")

# Load the neurodegenerative associations and get genes
gene_associations <- read.csv("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/Horvath_gene_enrichment_overlap_greaterThan1_neurodegenerative_associations.csv")

neuro_genes <- unique(gene_associations$gene)
cat(sprintf("Neurodegenerative genes: %s\n\n", paste(neuro_genes, collapse = ", ")))

# Load the hypermethylated CpGs
hyper_cpgs <- readRDS("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251208_hypermethylated_Horvath_CpGs_African_Ancestry.rds")

cat(sprintf("Total hypermethylated CpGs: %d\n", length(hyper_cpgs)))

# Build gene databases to get CpG-to-gene mappings
cat("\nBuilding gene databases...\n")
gene_dbs <- buildGeneDBs(hyper_cpgs, platform = "MSA")

cat(sprintf("Total genes in database: %d\n", length(gene_dbs)))

# Extract gene names from attributes
gene_names <- sapply(gene_dbs, function(x) attr(x, "gene_name"))

# Find which neurodegenerative genes are present
found_genes <- neuro_genes[neuro_genes %in% gene_names]

cat(sprintf("\nFound %d neurodegenerative genes in gene_dbs: %s\n", 
            length(found_genes), 
            paste(found_genes, collapse = ", ")))

# Get the CpGs for each found gene
if (length(found_genes) > 0) {
  cpg_gene_mapping <- data.frame()
  
  for (gene in found_genes) {
    # Find which element has this gene name
    idx <- which(gene_names == gene)
    
    for (i in idx) {
      cpgs <- gene_dbs[[i]]
      cat(sprintf("\n%s has %d CpGs\n", gene, length(cpgs)))
      
      cpg_gene_mapping <- rbind(cpg_gene_mapping,
                                data.frame(CpG = cpgs,
                                           gene = gene,
                                           stringsAsFactors = FALSE))
    }
  }
  
  # Aggregate by CpG (in case a CpG maps to multiple genes)
  neuro_cpg_mapping <- cpg_gene_mapping %>%
    group_by(CpG) %>%
    summarize(neuro_genes = paste(unique(gene), collapse = ";"), .groups = "drop")
  
  cat(sprintf("\nTotal unique CpGs mapping to neuro genes: %d\n", nrow(neuro_cpg_mapping)))
  
  # Find intersection with hypermethylated CpGs
  neuro_cpgs <- neuro_cpg_mapping$CpG
  intersection <- intersect(neuro_cpgs, hyper_cpgs)
  
  cat(sprintf("\nCpGs in neuro_cpg_mapping: %d\n", length(neuro_cpgs)))
  cat(sprintf("CpGs in hyper_cpgs: %d\n", length(hyper_cpgs)))
  cat(sprintf("Intersection: %d\n", length(intersection)))
  
  if (length(intersection) > 0) {
    # Filter to only the intersection
    neuro_hyper_cpgs <- neuro_cpg_mapping %>%
      filter(CpG %in% intersection)
    
    cat("\nCpGs that are both neurodegenerative-associated AND hypermethylated (with suffix):\n")
    print(neuro_hyper_cpgs)
    
    # Strip the suffix from CpG IDs to match the results file
    neuro_hyper_cpgs <- neuro_hyper_cpgs %>%
      mutate(CpG_base = gsub("_.*$", "", CpG))  # Remove everything after underscore
    
    cat("\nCpGs with base IDs (matching results format):\n")
    print(neuro_hyper_cpgs)
    
  } else {
    cat("\nNo overlap between neurodegenerative gene CpGs and hypermethylated CpGs\n")
    neuro_hyper_cpgs <- data.frame(CpG = character(0), neuro_genes = character(0), CpG_base = character(0))
  }
  
} else {
  cat("\nNone of the neurodegenerative genes found in gene_dbs\n")
  neuro_hyper_cpgs <- data.frame(CpG = character(0), neuro_genes = character(0), CpG_base = character(0))
}

# Annotate the results using CpG_base to match
if (nrow(neuro_hyper_cpgs) > 0) {
  results_annotated <- results %>%
    left_join(neuro_hyper_cpgs %>% select(CpG_base, neuro_genes), 
              by = c("CpG" = "CpG_base")) %>%
    mutate(
      has_neuro_gene = !is.na(neuro_genes)
    )
} else {
  results_annotated <- results %>%
    mutate(
      neuro_genes = NA_character_,
      has_neuro_gene = FALSE
    )
}

# Count CpGs
cat("\n=== CpG counts ===\n")
cat(sprintf("Total CpGs: %d\n", nrow(results_annotated)))
cat(sprintf("Neurodegenerative-associated CpGs: %d\n", sum(results_annotated$has_neuro_gene)))
print(table(results_annotated$Ancestry))

# Prepare labels - with both CpG and gene name
label_data <- results_annotated %>%
  filter(has_neuro_gene) %>%
  mutate(label = paste0(CpG, " (", neuro_genes, ")"))

cat(sprintf("\n=== Labeling %d CpGs ===\n", nrow(label_data)))
if (nrow(label_data) > 0) {
  cat("\nCpGs to be labeled:\n")
  print(label_data %>% select(CpG, neuro_genes, label, Ancestry, adj.P.Val, delta_beta))
}

# Define color palette - just use original Ancestry
ancestry_colors <- c(
  "Higher in AFR" = "#E69F00",  # orange
  "Higher in EUR" = "#009E73",  # bluish green
  "Not Significant" = "grey"
)

# Calculate plot limits
max_delta <- ceiling(max(abs(results_annotated$delta_beta * 100), na.rm = TRUE) / 5) * 5
legend_breaks <- c(0, 0.15, 0.3, 0.45, 0.6)

# Create volcano plot
volcano_neuro_annotated <- ggplot(results_annotated, 
                                  aes(x = delta_beta * 100, 
                                      y = -log10(adj.P.Val), 
                                      color = Ancestry, 
                                      size = abs_weight)) +
  geom_point(alpha = 0.6, shape = 16)

# Add highlighting and labels only if we have data
if (nrow(label_data) > 0) {
  volcano_neuro_annotated <- volcano_neuro_annotated +
    geom_point(data = label_data,
               aes(x = delta_beta * 100, y = -log10(adj.P.Val)),
               shape = 21, stroke = 1.5, fill = NA, color = "black", size = 5) +
    geom_text_repel(data = label_data,
                    aes(label = label),
                    size = 4.5,  # Increased from 3.5
                    color = "black",
                    fontface = "bold",
                    max.overlaps = 30,
                    box.padding = 1.2,  # Increased from 0.8
                    point.padding = 0.8,  # Increased from 0.5
                    segment.color = "black",
                    segment.size = 0.4,
                    min.segment.length = 0,
                    force = 5,  # Increased from 3
                    force_pull = 0.1,
                    nudge_x = 8,  # Increased from 5
                    direction = "y")
}

# Complete the plot
volcano_neuro_annotated <- volcano_neuro_annotated +
  scale_size_continuous(
    range = c(2, 6),
    breaks = legend_breaks,
    labels = function(x) sprintf("%.2f", x),
    limits = c(0, NA),
    name = "Clock Weight (absolute)"
  ) +
  scale_color_manual(values = ancestry_colors, name = "Methylation Level") +
  guides(
    color = guide_legend(override.aes = list(size = 3, alpha = 1), order = 1),
    size = guide_legend(override.aes = list(color = "black", alpha = 1), order = 2)
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  coord_cartesian(clip = "off") +
  annotate("text", x = -Inf, y = -log10(0.05), label = "FDR = 0.05",
           hjust = 1.05, vjust = 0.5, size = 14 / .pt, color = "black") +
  scale_x_continuous(breaks = seq(-max_delta, max_delta, by = 5),
                     limits = c(-max_delta-2, max_delta+2),
                     expand = expansion(mult = 0.02)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
  labs(x = expression(paste(Delta, "β(%)")),
       y = expression(-log[10](FDR)),
       title = "Differential Methylation with Neurodegenerative Gene Annotations") +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = unit(0.2, "cm"),
    plot.margin = margin(5.5, 20, 5.5, 20, "pt"),
    plot.title = element_text(size = 16, hjust = 0.5)
  )

print(volcano_neuro_annotated)

# Save
ggsave("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251209_Horvath_volcano_neurodegenerative_annotated.png", 
       plot = volcano_neuro_annotated, width = 12, height = 8, dpi = 300)

cat(sprintf("Hypermethylated CpGs mapping to neurodegenerative genes: %d\n", sum(results_annotated$has_neuro_gene, na.rm = TRUE)))