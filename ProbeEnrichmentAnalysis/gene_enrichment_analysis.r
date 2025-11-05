# This script runs a gene enrichment analysis, using KnowYourCG, for CpGs that are differently methylated by genetic ancestry
# Author: Rory Boyle rorytboyle@gmail.com
# Date: 05/11/2025

# Load required libraries
library(sesame)
library(sesameData)
library(knowYourCG)
library(dplyr)
library(ggplot2)

# Cache databases
sesameDataCache()

# Load your differential methylation results
dmp_cgs <- read.csv('/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251105_171_DiffMethylAnalysis_results.csv')

# Get significant CpGs (hypermethylated in African ancestry)
dmp_cpgs_sig <- dmp_cgs %>%
  filter(adj.P.Val < 0.05) %>%
  filter(Ancestry == "Higher in AFR")

query_cpgs <- dmp_cpgs_sig$CpG

# Map CpGs to genes and test enrichment
gene_dbs <- buildGeneDBs(query_cpgs, platform = "EPIC")
gene_res <- testEnrichment(query_cpgs, gene_dbs, platform = "EPIC")

# Get significant genes
sig_genes <- gene_res %>%
  as.data.frame() %>%
  filter(FDR < 0.05) %>%
  arrange(FDR)

cat(sprintf("Found %d significantly enriched genes (FDR < 0.05)\n", nrow(sig_genes)))

sig_genes_multi <- sig_genes %>%
  filter(overlap > 1) %>%
  arrange(p.value) %>%
  mutate(
    gene_label = paste0("Gene | ", gene_name),
    n_label = paste0("N = ", overlap),
    log2_OR = log2(estimate),
    log10_p = -log10(p.value),
    gene_order = factor(gene_label, levels = rev(gene_label))
  )

cat(sprintf("Genes with >1 CpG overlap: %d\n", nrow(sig_genes_multi)))

# Left panel: P-values with N labels embedded in bars
p_gene_left <- ggplot(sig_genes_multi, aes(x = gene_order, y = log10_p)) +
  geom_col(fill = "grey40", color = "black", width = 0.85, linewidth = 0.2) +
  # N labels inside bars (white text)
  geom_text(aes(label = n_label), 
            # hjust = 1.1 makes the label align right inside the bar (since the bar extends right)
            hjust = 1.1,
            size = 3,
            color = "white",
            fontface = "bold") +
  # Use a regular scale for -log10(P-value)
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)), # Ensures bars start at 0
    # Set breaks to correspond to the original P-value
    breaks = c(0, 2, 4, 6, 8),
    labels = c("1", expression(10^-2), expression(10^-4), expression(10^-6), expression(10^-8))
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = expression(-Log[10](italic(P)~value)) # Label the axis correctly
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 10),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    plot.margin = margin(5, 1, 5, 5)
  )

# Right panel: Log2(OR)
p_gene_right <- ggplot(sig_genes_multi, aes(x = gene_order, y = log2_OR)) +
  geom_col(fill = "grey40", color = "black", width = 0.85, linewidth = 0.2) +
  scale_x_discrete(position = "top") +
  coord_flip() +
  labs(
    x = NULL,
    y = expression(Log[2](OR))
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 10),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    plot.margin = margin(5, 5, 5, 1)
  )

# Combine gene panels
combined_genes <- p_gene_left + plot_spacer() + p_gene_right + 
  plot_layout(widths = c(1.3, 0.05, 1))

print(combined_genes)

# Export gene results
write.csv(sig_genes, 
          "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251105_gene_enrichment.csv",
          row.names = FALSE)