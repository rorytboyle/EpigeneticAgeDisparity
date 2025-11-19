# This script runs probe enrichment analyses for CpGs from a differential methylation analysis.
# Analyses include testing enrichment, using KnowYourCG, for chromatin states, transcription factor binding sites,
# CpG island context, and genomic regions, and testing genomic proximity of the CpGs.
# Author: Rory Boyle & Nadia Dehghani rorytboyle@gmail.com
# Date: 19th November 2025

library(sesame)
library(sesameData)
library(knowYourCG)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(gwasrapidd)

# Load differentially methylated cpgs ####
# Cache databases (only needed once)
sesameDataCache()

# Load your differential methylation results
dmp_cgs <- read.csv('/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251105_171_DiffMethylAnalysis_results.csv')

# Get significant CpGs (hypermethylated in African ancestry)
dmp_cpgs_sig <- dmp_cgs %>%
  filter(adj.P.Val < 0.05) %>%
  filter(Ancestry == "Higher in AFR")

query_cpgs <- dmp_cpgs_sig$CpG

cat(sprintf("Analyzing %d significant CpGs\n", length(query_cpgs)))

# Probe enrichment for chromatin states ####
chrom_res <- testEnrichment(query_cpgs, platform="EPIC", databases = "KYCG.EPIC.chromHMM.20211020")

# Chromatin state labels taken from here: https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html
# Create lookup table for chromatin states
chromhmm_labels <- data.frame(
  mnemonic = c("1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", 
               "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het", "10_TssBiv",
               "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies"),
  description = c("Active TSS", "Flanking Active TSS", "Transcr. at gene 5' and 3'",
                  "Strong transcription", "Weak transcription", "Genic enhancers",
                  "Enhancers", "ZNF genes & repeats", "Heterochromatin",
                  "Bivalent/Poised TSS", "Flanking Bivalent TSS/Enh",
                  "Bivalent Enhancer", "Repressed PolyComb",
                  "Weak Repressed PolyComb", "Quiescent/Low"),
  color_name = c("Red", "Orange Red", "LimeGreen", "Green", "DarkGreen",
                 "GreenYellow", "Yellow", "Medium Aquamarine", "PaleTurquoise",
                 "IndianRed", "DarkSalmon", "DarkKhaki", "Silver", "Gainsboro", "White"),
  color_code = c("255,0,0", "255,69,0", "50,205,50", "0,128,0", "0,100,0",
                 "194,225,5", "255,255,0", "102,205,170", "138,145,208",
                 "205,92,92", "233,150,122", "189,183,107", "128,128,128",
                 "192,192,192", "255,255,255"),
  stringsAsFactors = FALSE
)

# Convert color codes to hex
chromhmm_labels$hex_color <- sapply(strsplit(chromhmm_labels$color_code, ","), function(x) {
  rgb(as.numeric(x[1]), as.numeric(x[2]), as.numeric(x[3]), maxColorValue = 255)
})

# Filter and prepare chromatin state data with clean labels
chrom_sig <- chrom_res %>%
  as.data.frame() %>%
  # filter(FDR < 0.05) %>%
  arrange(p.value) %>%
  left_join(chromhmm_labels, by = c("dbname" = "mnemonic")) %>%
  mutate(
    state_num = sub("_.*", "", dbname),
    state_abbrev = sub(".*_", "", dbname),
    clean_label = paste0(description, " (", state_num, "_", state_abbrev, ")"),
    log2_OR = log2(estimate),
    log10_p = -log10(p.value),
    state_order = factor(clean_label, levels = rev(clean_label)),
    font_face = ifelse(FDR < 0.05, "bold", "plain")
  )

# Two-panel plot with chromHMM colors
# Two-panel plot with chromHMM colors
p_chrom_left <- ggplot(chrom_sig, aes(x = state_order, y = log10_p)) +
  geom_col(aes(fill = hex_color), color = "black", width = 0.85, linewidth = 0.2) +
  scale_fill_identity() +
  geom_text(aes(label = paste0("N = ", overlap)), 
            hjust = 1.1, size = 3, color = "black", fontface = "bold") +
  
  # Vertical line for P-value < 0.05
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.5) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    breaks = seq(0, ceiling(max(chrom_sig$log10_p)), by = 2)
  ) +
  coord_flip() +
  labs(x = NULL, y = expression(-Log[10](italic(P)~value))) +
  theme_classic(base_size = 11) +
  theme(
    # Use the font_face vector for the label order
    axis.text.y = element_text(size = 9, color = "black", face = chrom_sig$font_face[match(levels(chrom_sig$state_order), chrom_sig$clean_label)]),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 10),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    plot.margin = margin(5, 1, 5, 5)
  )

p_chrom_right <- ggplot(chrom_sig, aes(x = state_order, y = log2_OR)) +
  geom_col(aes(fill = hex_color), color = "black", width = 0.85, linewidth = 0.2) +
  scale_fill_identity() +
  coord_flip() +
  labs(x = NULL, y = expression(Log[2](OR))) +
  theme_classic(base_size = 11) +
  theme(
    # also apply the bolding to the right plot's y-axis labels  but since they are blanked out, 
    # it won't show.
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

combined_chrom <- p_chrom_left + plot_spacer() + p_chrom_right + 
  plot_layout(widths = c(1.3, 0.05, 1))

print(combined_chrom)

ggsave("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/chromatin_state_enrichment.png", combined_chrom, width = 10, height = 6, dpi = 300)

# Probe enrichment for transcription factor binding sites ####
tfbs_res <- testEnrichment(query_cpgs, platform="EPIC", databases = "KYCG.EPIC.TFBSconsensus.20211013") 

tfbs_sig <- tfbs_res %>%
  as.data.frame() %>%
  filter(p.value < 0.05) %>% # No sig results at FDR so just show nominally sig
  arrange(p.value) %>%
  # slice_head(n = 25) %>% # Top 25 for visibility
  mutate(
    tf_label = dbname, 
    n_label = paste0("N = ", overlap),
    log2_OR = log2(estimate),
    log10_p = -log10(p.value),
    tf_order = factor(tf_label, levels = rev(tf_label))
  )

# Left panel: P-values with N labels embedded in bars
p_tfbs_left <- ggplot(tfbs_sig, aes(x = tf_order, y = log10_p)) +
  geom_col(fill = "grey40", color = "black", width = 0.85, linewidth = 0.2) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed", linewidth = 0.8) + # could remove for consistency
  geom_text(aes(label = n_label), 
            hjust = 1.1,
            size = 3,
            color = "white",
            fontface = "bold") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    breaks = c(0, 2, 4, 6, 8, 10),
    labels = c("1", expression(10^-2), expression(10^-4), expression(10^-6), 
               expression(10^-8), expression(10^-10))
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = expression(-Log[10](italic(P)~value))
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
p_tfbs_right <- ggplot(tfbs_sig, aes(x = tf_order, y = log2_OR)) +
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

# Combine TFBS panels
combined_tfbs <- p_tfbs_left + plot_spacer() + p_tfbs_right + 
  plot_layout(widths = c(1.3, 0.05, 1)) +
  plot_annotation(
    title = "Probe enrichment for transcription factor binding sites"
  )

print(combined_tfbs)

ggsave("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/transcription_factor_binding_site_enrichment.png", combined_tfbs, width = 10, height = 6, dpi = 300)

# No enrichment for transcription factor binding sites at FDR (only nominally significant)

# Probe enrichment for CpG context ####
cgi_res <- testEnrichment(query_cpgs, platform="EPIC", databases = "KYCG.EPIC.CGI.20210713")

cgi_sig <- cgi_res %>%
  as.data.frame() %>%
  # filter(FDR < 0.05) %>% # No sig results at FDR so just show all
  mutate(
    category = case_when(
      grepl("island", dbname, ignore.case = TRUE) ~ "Island",
      grepl("shore", dbname, ignore.case = TRUE) ~ "Shore", 
      grepl("shelf", dbname, ignore.case = TRUE) ~ "Shelf",
      # grepl("open_sea", dbname, ignore.case = TRUE) ~ "Open Sea", #not actually in results
      TRUE ~ dbname
    ),
    log10_p = -log10(p.value),
    log2_OR = log2(estimate)
  )

# Prepare CGI enrichment data with proper labels
cgi_annotated <- cgi_res %>%
  as.data.frame() %>%
  mutate(
    context = case_when(
      grepl("N_Shelf", dbname, ignore.case = TRUE) ~ "N. Shelf",
      grepl("N_Shore", dbname, ignore.case = TRUE) ~ "N. Shore",
      grepl("island", dbname, ignore.case = TRUE) & !grepl("shore|shelf", dbname, ignore.case = TRUE) ~ "CGI",
      grepl("S_Shore", dbname, ignore.case = TRUE) ~ "S. Shore",
      grepl("S_Shelf", dbname, ignore.case = TRUE) ~ "S. Shelf",
      # grepl("sea", dbname, ignore.case = TRUE) ~ "Open Sea", # not in results
      TRUE ~ dbname
    ),
    log2_OR = log2(estimate),
    log10_p = -log10(p.value)
  )

# Create schematic layout
schematic_layout <- data.frame(
  context = factor(c("N. Shelf", "N. Shore", "CGI", "S. Shore", "S. Shelf"),
                   levels = c("N. Shelf", "N. Shore", "CGI", "S. Shore", "S. Shelf")),
  x_start = c(0, 2, 4, 6, 8),
  x_end = c(2, 4, 6, 8, 10),
  y_bottom = 0,
  y_top = 1,
  # Colors from your screenshot (blue to green gradient with yellow/orange center)
  region_color = c("#4575b4", "#91bfdb", "#ffffbf", "#fc8d59", "#d73027")
)

# Merge with enrichment data
plot_data <- schematic_layout %>%
  left_join(cgi_annotated, by = "context") %>%
  mutate(
    x_mid = (x_start + x_end) / 2,
    significant = !is.na(FDR) & FDR < 0.05
  )

# Create the plot
cpg_context_plot <- ggplot(plot_data) +
  # Draw rectangles with colors from screenshot
  geom_rect(aes(xmin = x_start, xmax = x_end, 
                ymin = y_bottom, ymax = y_top, 
                fill = region_color),
            color = "black", linewidth = 1.2) +
  scale_fill_identity() +  # Use the exact colors specified
  # Add significance borders for FDR < 0.05
  geom_rect(data = filter(plot_data, significant),
            aes(xmin = x_start, xmax = x_end, 
                ymin = y_bottom, ymax = y_top),
            color = "red", fill = NA, linewidth = 2.5) +
  # Add context labels at top
  geom_text(aes(x = x_mid, y = 1.15, label = context), 
            size = 5, fontface = "bold") +
  # Add N (overlap)
  geom_text(aes(x = x_mid, y = 0.75, 
                label = sprintf("N = %d", ifelse(is.na(overlap), 0, overlap))), 
            size = 4, fontface = "bold") +
  # Add Log2(OR) - always shown with 3 decimal places
  geom_text(aes(x = x_mid, y = 0.55, 
                label = sprintf("Log2(OR) = %.3f", 
                                ifelse(is.na(log2_OR), 0, log2_OR))), 
            size = 3.5) +
  # Add FDR - always shown with 3 decimal places, with color coding
  geom_text(aes(x = x_mid, y = 0.35, 
                label = sprintf("FDR = %.3f", ifelse(is.na(FDR), 1, FDR)),
                color = significant), 
            size = 3.5, fontface = "bold") +
  # Add p-value - always shown with 3 decimal places
  geom_text(aes(x = x_mid, y = 0.15, 
                label = sprintf("p = %.3f", ifelse(is.na(p.value), 1, p.value))), 
            size = 3) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "black"),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(-0.5, 10.5), ylim = c(0, 1.3), clip = "off") +
  labs(title = "CpG Island Context Enrichment Analysis") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.margin = margin(20, 20, 20, 20)
  )

cpg_context_plot
# Save high-res version
ggsave("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/cpg_island_context_enrichment.png", cpg_context_plot, width = 10, height = 6, dpi = 300)

# Probe genomic proximity of CpGs ####
# https://github.com/zhou-lab/knowYourCG/blob/devel/R/proximal_genes.R # This won't run
proximity_res <- testProbeProximity(query_cpgs, platform="EPIC")
proximity_res

# Visualize

# Step 1 Get annotations for CpGs
anno <- sesameData_getManifestGRanges("EPIC")
query_anno <- anno[names(anno) %in% query_cpgs]

# Convert to data frame
cpg_positions <- as.data.frame(query_anno) %>%
  mutate(CpG = names(query_anno)) %>%
  select(CpG, seqnames, start)

# Identify cluster members
cluster_cpgs <- cpg_positions %>%
  filter(seqnames == "chr20", 
         start >= 43945701, 
         start <= 43945985)

cat("CpGs in the cluster:\n")
print(cluster_cpgs)

# Add cluster membership to all CpGs
cpg_positions <- cpg_positions %>%
  mutate(in_cluster = CpG %in% cluster_cpgs$CpG)

# Step 2 Plot CpG distribution across chromosomes
# Create complete chromosome list
all_chroms <- paste0("chr", c(1:22, "X", "Y"))

cpg_summary <- cpg_positions %>%
  group_by(seqnames) %>%
  dplyr::summarize(
    n_cpgs = dplyr::n(),
    has_cluster = any(in_cluster),
    .groups = "drop"
  ) %>%
  # Add missing chromosomes with 0 CpGs
  complete(seqnames = all_chroms, fill = list(n_cpgs = 0, has_cluster = FALSE)) %>%
  mutate(seqnames = factor(seqnames, levels = all_chroms))

# Manually set chr_x to be chromosome with identified cluster from proximity_res
chr_x <- 20

# Prepare cluster annotation data with specific y positions
cluster_annotation <- cluster_cpgs %>%
  mutate(
    seqnames = factor(seqnames, levels = all_chroms),
    label = CpG
  ) %>%
  left_join(cpg_summary %>% select(seqnames, n_cpgs), by = "seqnames") %>%
  arrange(start) %>%
  mutate(
    label_y = seq(7.5, 15, length.out = dplyr::n())
  )

# Calculate distance annotation
if (nrow(cluster_annotation) == 2) {
  distance_bp <- abs(diff(cluster_annotation$start))
  distance_y <- mean(cluster_annotation$label_y)
  distance_label <- paste0(distance_bp, " bp")
}

cpg_genomic_proximity <- ggplot(cpg_summary, aes(x = seqnames, y = n_cpgs, fill = has_cluster)) +
  geom_col() +
  # Add line segments from bar to labels
  geom_segment(data = cluster_annotation,
               aes(x = seqnames, xend = seqnames, 
                   y = n_cpgs, yend = label_y - 0.5),
               inherit.aes = FALSE,
               color = "red", linewidth = 0.5, linetype = "dashed") +
  # Add text labels for CpGs (centered)
  geom_text(data = cluster_annotation,
            aes(x = seqnames, y = label_y, label = label),
            inherit.aes = FALSE,
            size = 3, fontface = "bold", color = "red") +
  # Add distance label (horizontal, to the left)
  annotate("text", 
           x = chr_x - 2, 
           y = distance_y,
           label = distance_label,
           size = 3.5, fontface = "italic", color = "darkred", hjust = 1) +
  # Add arrow to UPPER CpG (stop BELOW the text)
  annotate("segment",
           x = chr_x - 1.9,
           xend = chr_x,
           y = distance_y,
           yend = cluster_annotation$label_y[2] - 0.3,
           color = "darkred", linewidth = 0.6, 
           arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  # Add arrow to LOWER CpG (stop ABOVE the text)
  annotate("segment",
           x = chr_x - 1.9,
           xend = chr_x,
           y = distance_y,
           yend = cluster_annotation$label_y[1] + 0.3,
           color = "darkred", linewidth = 0.6,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  scale_fill_manual(values = c("FALSE" = "gray60", "TRUE" = "red"),
                    name = "Contains cluster") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Chromosome", y = "Number of CpGs",
       title = "Distribution of Hypermethylated CpGs Across Chromosomes",
       subtitle = paste0("Spatial clustering P = ", format(proximity_res$Stats$P.val, scientific = TRUE))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(cpg_genomic_proximity)

ggsave("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/cpg_genomic_proximity.png", cpg_genomic_proximity, width = 10, height = 6, dpi = 300)


# Step 3 Get genes near the cluster CpGs
if (nrow(cluster_cpgs) > 0) {
  gene_dbs <- buildGeneDBs(cluster_cpgs$CpG, platform = "EPIC")
  gene_res <- testEnrichment(cluster_cpgs$CpG, gene_dbs, platform = "EPIC")
  
  cluster_genes <- gene_res %>%
    as.data.frame() %>%
    filter(overlap > 0) %>%
    arrange(desc(overlap))
  
  cat("\nGenes associated with cluster CpGs:\n")
  print(cluster_genes %>% select(gene_name, overlap, estimate, p.value))
}


# Probe enrichment for genomic region ####
# no significant enrichment found
metagene_res <- testEnrichment(
  query_cpgs,
  platform = "EPIC",
  databases = "KYCG.EPIC.metagene.20220126"
)
metagene_df <- metagene_res %>%
  as.data.frame() %>%
  arrange(FDR)
cat("Gene region enrichment:\n")
print(metagene_df %>% select(dbname, estimate, overlap, p.value, FDR))

# Get the actual labels from the database
metagene_db <- getDBs("EPIC.metagene")
metagene_actual_labels <- sapply(metagene_db, function(x) attr(x, "label"))

# Add to your dataframe
metagene_df <- metagene_res %>%
  as.data.frame() %>%
  arrange(FDR) %>%
  mutate(region_label = metagene_actual_labels[as.character(dbname)])

# Visualize gene regions with actual labels
p_metagene <- ggplot(metagene_df, aes(x = reorder(region_label, as.numeric(dbname)), 
                                      y = -log10(FDR))) +
  geom_col(fill = "gray60", color = "black", width = 0.7) +
  geom_text(aes(label = overlap), vjust = -0.5, size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(
    title = "Gene Region Enrichment",
    # subtitle = "Location relative to gene structure",
    x = "Gene Region",
    y = "-log10(FDR)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

print(p_metagene)
ggsave("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/enrichment_for_genomic_region.png", p_metagene, width = 10, height = 6, dpi = 300)

