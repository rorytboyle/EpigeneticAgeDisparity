# This script runs probe enrichment analyses for CpGs from a differential methylation analysis on Horvath clock CpGs by genetic ancestry.
# Analyses include testing enrichment, using KnowYourCG, for chromatin states, transcription factor binding sites, tissue signature,
# CpG island context, and genomic regions, and testing genomic proximity of the CpGs.
# Author: Rory Boyle & Nadia Dehghani rorytboyle@gmail.com
# Date: 19th November 2025
# Updated: 1st December 2025 to use development version of KnowYourCG and use MSA platform
# Updated: 8th December 2025 to fix text label positioning for narrow bars with dynamic threshold & fix bug in transcription factor binding site enrichment analysis

# Set up ####
# Note: This requires the development version of knowYourCG from GitHub: BiocManager::install('zhou-lab/knowYourCG')
# This may require clean up before you can install successfully. If you get an error installing the dev version, try below:
# # 1. Find all potentially corrupted packages
# lib_path <- .libPaths()[1]
# list.files(lib_path)
# # 2. Remove the known corrupted packages
# remove.packages(c("BiocManager", "sesameData", "knowYourCG"))
# # 3. Reinstall from scratch
# install.packages("BiocManager")
# BiocManager::install("sesameData")
# BiocManager::install('zhou-lab/knowYourCG')
# # 4. Load library
# library(knowYourCG)

library(sesame)
library(sesameData)
library(knowYourCG)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(gwasrapidd)

# Set seed for reproducibility
set.seed(123)

# Load differentially methylated cpgs ####
# Cache databases (only needed once)
sesameDataCache()

# Load your differential methylation results 
# Read in CpGs that are significantly hypermethylated in African Ancestry at FDR < 0.05 when adjusted for cell type proportion analyses
# ~/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/code/prep_CpG_query_Horvath.R
query_cpgs <- readRDS("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251208_hypermethylated_Horvath_CpGs_African_Ancestry.rds")

cat(sprintf("Analyzing %d significant CpGs\n", length(query_cpgs)))

# Create helper function to ensure that overlap labels always correctly plotted on bar plots
# Dynamic threshold calculation function
calculate_label_threshold <- function(log10_p_values, min_bar_width_fraction = 0.15) {
  # Calculate threshold as a fraction of the maximum -log10(p) value
  # This ensures text fits inside bars that are at least X% of the max height
  max_log10p <- max(log10_p_values, na.rm = TRUE)
  threshold <- max_log10p * min_bar_width_fraction
  return(threshold)
}


# Probe enrichment for chromatin states ####
chrom_res <- testEnrichment(query_cpgs, platform="MSA", databases = "MSA.ChromHMM.20220303")

# Annotate chromHMM results with state names from development version of KnowYourCG
# in previous versions, dbname is returned in chrom_res but in the dev version using MSA, it is not. 
# So these states have to be annotated manually based on the known order of states in the MSA ChromHMM database.
annotate_chromhmm_results <- function(enrichment_results, 
                                      platform = "MSA",
                                      db_path = "~/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/MSA_KYCG_Knowledgebases/KYCG.MSA.ChromHMM.20220303.rda") {
  
  msa_chromhmm_states <- data.frame(
    chromatin_state = as.character(1:19),
    state_name = c(
      "EnhA1", "EnhA2", "EnhBiv", "EnhG1", "EnhG2", "EnhWk", "Het", "NA", "Quies",
      "ReprPC", "ReprPCWk", "TssA", "TssBiv", "TssFlnk", "TssFlnkD", "TssFlnkU",
      "Tx", "TxWk", "ZNF/Rpts"
    ),
    description = c(
      "Active Enhancer 1",
      "Active Enhancer 2", 
      "Bivalent Enhancer",
      "Genic enhancer 1",
      "Genic enhancer 2",
      "Weak Enhancer",
      "Heterochromatin",
      "Not Annotated",
      "Quiescent/Low",
      "Repressed PolyComb",
      "Weak Repressed PolyComb",
      "Active TSS",
      "Bivalent/Poised TSS",
      "Flanking TSS",
      "Flanking TSS Downstream",
      "Flanking TSS Upstream",
      "Strong transcription",
      "Weak transcription",
      "ZNF genes & repeats"
    ),
    stringsAsFactors = FALSE
  )
  
  # Add chromatin state from rownames
  enrichment_results$chromatin_state <- rownames(enrichment_results)
  
  # Merge with state names
  annotated <- merge(enrichment_results, msa_chromhmm_states, 
                     by = "chromatin_state", all.x = TRUE)
  
  # Sort by FDR
  annotated <- annotated[order(annotated$FDR), ]
  
  # Reorder columns for better readability
  col_order <- c("chromatin_state", "state_name", "description", "overlap", "nD", "nQ",
                 "estimate", "p.value", "FDR", "cf_overlap", "cf_Jaccard")
  col_order <- col_order[col_order %in% names(annotated)]
  other_cols <- setdiff(names(annotated), col_order)
  annotated <- annotated[, c(col_order, other_cols)]
  
  return(annotated)
}

# Re-run with corrected annotations
chrom_res_annotated <- annotate_chromhmm_results(chrom_res)

# MSA ChromHMM state labels using expanded 18-state Roadmap model colors
# Reference: https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html
# Note: MSA database has its own numbering (1-19) that doesn't match Roadmap
msa_chromhmm_labels <- data.frame(
  state_num = as.character(1:19),
  state_abbrev = c("EnhA1", "EnhA2", "EnhBiv", "EnhG1", "EnhG2", "EnhWk", 
                   "Het", "NA", "Quies", "ReprPC", "ReprPCWk", "TssA", "TssBiv", 
                   "TssFlnk", "TssFlnkD", "TssFlnkU", "Tx", "TxWk", "ZNF/Rpts"),
  description = c("Active Enhancer 1", "Active Enhancer 2", "Bivalent Enhancer",
                  "Genic enhancer 1", "Genic enhancer 2", "Weak Enhancer",
                  "Heterochromatin", "Not Annotated", "Quiescent/Low",
                  "Repressed PolyComb", "Weak Repressed PolyComb", "Active TSS",
                  "Bivalent/Poised TSS", "Flanking TSS", "Flanking TSS Downstream",
                  "Flanking TSS Upstream", "Strong transcription", 
                  "Weak transcription", "ZNF genes & repeats"),
  # Roadmap 18-state model colors (mapped by function, not by number)
  color_name = c("Orange", "Orange", "DarkKhaki", "GreenYellow", "GreenYellow",
                 "Yellow", "PaleTurquoise", "Gray", "White", "Silver",
                 "Gainsboro", "Red", "IndianRed", "OrangeRed", "OrangeRed",
                 "OrangeRed", "Green", "DarkGreen", "MediumAquamarine"),
  color_code = c("255,195,77", "255,195,77", "189,183,107", "194,225,5", "194,225,5",
                 "255,255,0", "138,145,208", "128,128,128", "255,255,255", "128,128,128",
                 "192,192,192", "255,0,0", "205,92,92", "255,69,0", "255,69,0",
                 "255,69,0", "0,128,0", "0,100,0", "102,205,170"),
  stringsAsFactors = FALSE
)

# Convert color codes to hex
msa_chromhmm_labels$hex_color <- sapply(strsplit(msa_chromhmm_labels$color_code, ","), function(x) {
  rgb(as.numeric(x[1]), as.numeric(x[2]), as.numeric(x[3]), maxColorValue = 255)
})

# Add chromatin state from rownames to results
chrom_res$chromatin_state <- rownames(chrom_res)

# Filter and prepare chromatin state data with clean labels (NO numeric prefix)
chrom_sig <- chrom_res %>%
  as.data.frame() %>%
  left_join(msa_chromhmm_labels, by = c("chromatin_state" = "state_num")) %>%
  arrange(p.value) %>%
  mutate(
    clean_label = paste0(description, " (", state_abbrev, ")"),
    log2_OR = estimate,
    log10_p = -log10(p.value),
    state_order = factor(clean_label, levels = rev(clean_label)),
    font_face = ifelse(FDR < 0.05, "bold", "plain"),
    log10_p_plot = ifelse(overlap > 0, log10_p, NA),
    log2_OR_plot = ifelse(overlap > 0, log2_OR, NA)
  )

# Two-panel plot with chromHMM colors
chrom_threshold <- calculate_label_threshold(chrom_sig$log10_p_plot)

p_chrom_left <- ggplot(chrom_sig, aes(x = state_order, y = log10_p_plot)) +
  geom_col(aes(fill = hex_color), color = "black", width = 0.85, linewidth = 0.2) +
  scale_fill_identity() +
  # Text inside bar (for wider bars)
  geom_text(data = filter(chrom_sig, overlap > 0, log10_p_plot >= chrom_threshold),
            aes(label = paste0("N = ", overlap)), 
            hjust = 1.1, size = 3, color = "black", fontface = "bold") +
  # Text outside bar (for narrow bars)
  geom_text(data = filter(chrom_sig, overlap > 0, log10_p_plot < chrom_threshold),
            aes(label = paste0("N = ", overlap)), 
            hjust = -0.1, size = 3, color = "black", fontface = "bold") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.5) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    breaks = seq(0, ceiling(max(chrom_sig$log10_p, na.rm = TRUE)), by = 0.5)
  ) +
  coord_flip() +
  labs(x = NULL, y = expression(-Log[10](italic(P)~value))) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 9, color = "black", 
                               face = chrom_sig$font_face[match(levels(chrom_sig$state_order), 
                                                                chrom_sig$clean_label)]),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 10),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    plot.margin = margin(5, 1, 5, 5)
  )

p_chrom_right <- ggplot(chrom_sig, aes(x = state_order, y = log2_OR_plot)) +
  geom_col(aes(fill = hex_color), color = "black", width = 0.85, linewidth = 0.2) +
  scale_fill_identity() +
  coord_flip() +
  labs(x = NULL, y = expression(Log[2](OR))) +
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

combined_chrom <- p_chrom_left + plot_spacer() + p_chrom_right + 
  plot_layout(widths = c(1.3, 0.05, 1))

print(combined_chrom)

ggsave("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/Horvath_chromatin_state_enrichment.png", combined_chrom, width = 10, height = 6, dpi = 300)

write.csv(chrom_res_annotated, 
          "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/Horvath_chromatin_state_enrichment_results.csv",
          row.names = FALSE)

# Probe enrichment for transcription factor binding sites ####
tfbs_res <- testEnrichment(query_cpgs, platform="MSA", databases = "KYCG.MSA.TFBSrm.20221005") 

# Annotate transcription factor binding sites
annotate_tfbs_results <- function(enrichment_results,
                                  db_path = "~/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/MSA_KYCG_Knowledgebases/KYCG.MSA.TFBSrm.20221005.rda") {
  
  # Load the database
  load(db_path, envir = environment())
  db_obj <- KYCG.MSA.TFBSrm.20221005
  
  # Extract TF names from dbname attributes
  tfbs_names <- sapply(db_obj, function(x) {
    gsub("TFBSrm;", "", attr(x, "dbname"))
  })
  
  # Create mapping dataframe
  tfbs_mapping <- data.frame(
    tfbs_index = as.character(1:length(tfbs_names)),
    tf_name = tfbs_names,
    stringsAsFactors = FALSE
  )
  
  # Add TFBS index from rownames
  enrichment_results$tfbs_index <- rownames(enrichment_results)
  
  # Merge with TF names
  annotated <- merge(enrichment_results, tfbs_mapping, 
                     by = "tfbs_index", all.x = TRUE)
  
  # Sort by p-value
  annotated <- annotated[order(annotated$p.value), ]
  
  return(annotated)
}

tfbs_annotated <- annotate_tfbs_results(tfbs_res)

tfbs_sig <- tfbs_annotated %>%
  as.data.frame() %>%
  filter(p.value < 0.05) %>% # plot nominally significant (bolded by FDR)
  filter(overlap > 1) %>%
  arrange(p.value) %>%
  slice_head(n = 25) %>% # Top 25 for visibility
  mutate(
    tf_label = tf_name, 
    n_label = paste0("N = ", overlap),
    log2_OR = log2(estimate),
    log10_p = -log10(p.value),
    tf_order = factor(tf_label, levels = rev(tf_label)),
    font_face = ifelse(FDR < 0.05, "bold", "plain")
  )

# Check if there are any significant results
if (nrow(tfbs_sig) > 0) {
  
  tfbs_threshold <- calculate_label_threshold(tfbs_sig$log10_p)
  
  p_tfbs_left <- ggplot(tfbs_sig, aes(x = tf_order, y = log10_p)) +
    geom_col(fill = "grey40", color = "black", width = 0.85, linewidth = 0.2) +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed", linewidth = 0.8) +
    # Text inside bar (for wider bars)
    geom_text(data = filter(tfbs_sig, log10_p >= tfbs_threshold),
              aes(label = n_label), 
              hjust = 1.1, size = 3, color = "white", fontface = "bold") +
    # Text outside bar (for narrow bars)
    geom_text(data = filter(tfbs_sig, log10_p < tfbs_threshold),
              aes(label = n_label), 
              hjust = -0.1, size = 3, color = "black", fontface = "bold") +
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
      axis.text.y = element_text(size = 9, color = "black",
                                 face = tfbs_sig$font_face[match(levels(tfbs_sig$tf_order), 
                                                                 tfbs_sig$tf_label)]),
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
      title = "Probe enrichment for transcription factor binding sites (top 25)"
    )
  
  print(combined_tfbs)
  
  ggsave("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/Horvath_transcription_factor_binding_site_enrichment.png", 
         combined_tfbs, width = 10, height = 6, dpi = 300)
  
} else {
  cat("No transcription factors with FDR < 0.05 and overlap > 1\n")
  cat("Top transcription factors by p-value:\n")
  print(tfbs_annotated %>% 
          filter(overlap > 0) %>% 
          arrange(p.value) %>% 
          select(tf_name, overlap, estimate, p.value, FDR) %>%
          head(10))
}

# Save all TFBS results regardless
write.csv(tfbs_annotated, 
          "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/Horvath_tfbs_enrichment_results.csv",
          row.names = FALSE)

# Probe enrichment for tissue signature (TiSigLoyfer) ####
tissue_res <- testEnrichment(query_cpgs, platform="MSA", databases = "KYCG.MSA.TiSigLoyfer.20221209")

# Load the tissue signature database to get tissue names
load("~/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/MSA_KYCG_Knowledgebases/KYCG.MSA.TiSigLoyfer.20221209.rda")

# Extract tissue names
tissue_names <- sapply(KYCG.MSA.TiSigLoyfer.20221209, function(x) {
  gsub("TiSigLoyfer;", "", attr(x, "dbname"))
})

# Create mapping dataframe
tissue_mapping <- data.frame(
  tissue_index = as.character(1:length(tissue_names)),
  tissue_name = tissue_names,
  stringsAsFactors = FALSE
)

# Add tissue index from rownames
tissue_res$tissue_index <- rownames(tissue_res)

# Annotate tissue results
tissue_annotated <- tissue_res %>%
  as.data.frame() %>%
  left_join(tissue_mapping, by = "tissue_index") %>%
  arrange(p.value)

# Filter for significant or nominally significant results
tissue_sig <- tissue_annotated %>%
  filter(p.value < 0.05) %>%
  filter(overlap > 1) %>%
  arrange(p.value) %>%
  mutate(
    tissue_label = tissue_name,
    n_label = paste0("N = ", overlap),
    log2_OR = estimate,
    log10_p = -log10(p.value),
    tissue_order = factor(tissue_label, levels = rev(tissue_label)),
    font_face = ifelse(FDR < 0.05, "bold", "plain")
  )

# Check if there are any significant results
if (nrow(tissue_sig) > 0) {
  
  tissue_threshold <- calculate_label_threshold(tissue_sig$log10_p)
  
  p_tissue_left <- ggplot(tissue_sig, aes(x = tissue_order, y = log10_p)) +
    geom_col(fill = "steelblue", color = "black", width = 0.85, linewidth = 0.2) +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed", linewidth = 0.8) +
    # Text inside bar (for wider bars)
    geom_text(data = filter(tissue_sig, log10_p >= tissue_threshold),
              aes(label = n_label), 
              hjust = 1.1, size = 3, color = "white", fontface = "bold") +
    # Text outside bar (for narrow bars)
    geom_text(data = filter(tissue_sig, log10_p < tissue_threshold),
              aes(label = n_label), 
              hjust = -0.1, size = 3, color = "black", fontface = "bold") +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05))
    ) +
    coord_flip() +
    labs(
      x = NULL,
      y = expression(-Log[10](italic(P)~value))
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 9, color = "black",
                                 face = tissue_sig$font_face[match(levels(tissue_sig$tissue_order), 
                                                                   tissue_sig$tissue_label)]),
      axis.text.x = element_text(size = 9, color = "black"),
      axis.title.x = element_text(size = 10),
      panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor.x = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      plot.margin = margin(5, 1, 5, 5)
    )
  
  # Right panel: Log2(OR)
  p_tissue_right <- ggplot(tissue_sig, aes(x = tissue_order, y = log2_OR)) +
    geom_col(fill = "steelblue", color = "black", width = 0.85, linewidth = 0.2) +
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
  
  # Combine tissue panels
  combined_tissue <- p_tissue_left + plot_spacer() + p_tissue_right + 
    plot_layout(widths = c(1.3, 0.05, 1)) +
    plot_annotation(
      title = "Tissue Signature Enrichment (TiSigLoyfer)"
    )
  
  print(combined_tissue)
  
  # Save plot
  ggsave("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/Horvath_tissue_signature_enrichment_TiSigLoyfer.png", 
         combined_tissue, width = 10, height = 6, dpi = 300)
  
} else {
  cat("No significant tissue enrichment (TiSigLoyfer) found at p < 0.05 with overlap > 1\n")
  cat("Top tissues by p-value:\n")
  print(tissue_annotated %>% 
          filter(overlap > 0) %>% 
          arrange(p.value) %>% 
          select(tissue_name, overlap, estimate, p.value, FDR) %>%
          head(10))
}

# Save tissue enrichment results
write.csv(tissue_annotated,
          "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/Horvath_tissue_signature_enrichment_TiSigLoyfer_results.csv",
          row.names = FALSE)

# Probe enrichment for tissue signature (TiSigBLUEPRINT) ####
tissue_res <- testEnrichment(query_cpgs, platform="MSA", databases = "KYCG.MSA.TiSigBLUEPRINT.20221209")

# Load the tissue signature database to get tissue names
load("~/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/MSA_KYCG_Knowledgebases/KYCG.MSA.TiSigBLUEPRINT.20221209.rda")

# Extract tissue names
tissue_names <- sapply(KYCG.MSA.TiSigBLUEPRINT.20221209, function(x) {
  gsub("TiSigBLUEPRINT;", "", attr(x, "dbname"))
})

# Create mapping dataframe
tissue_mapping <- data.frame(
  tissue_index = as.character(1:length(tissue_names)),
  tissue_name = tissue_names,
  stringsAsFactors = FALSE
)

# Add tissue index from rownames
tissue_res$tissue_index <- rownames(tissue_res)

# Annotate tissue results
tissue_annotated <- tissue_res %>%
  as.data.frame() %>%
  left_join(tissue_mapping, by = "tissue_index") %>%
  arrange(p.value)

# Filter for significant or nominally significant results
tissue_sig <- tissue_annotated %>%
  filter(p.value < 0.05) %>%
  filter(overlap > 1) %>%
  arrange(p.value) %>%
  mutate(
    tissue_label = tissue_name,
    n_label = paste0("N = ", overlap),
    log2_OR = estimate,
    log10_p = -log10(p.value),
    tissue_order = factor(tissue_label, levels = rev(tissue_label)),
    font_face = ifelse(FDR < 0.05, "bold", "plain")
  )

# Check if there are any significant results
if (nrow(tissue_sig) > 0) {
  
  tissue_threshold <- calculate_label_threshold(tissue_sig$log10_p)
  
  p_tissue_left <- ggplot(tissue_sig, aes(x = tissue_order, y = log10_p)) +
    geom_col(fill = "steelblue", color = "black", width = 0.85, linewidth = 0.2) +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed", linewidth = 0.8) +
    # Text inside bar (for wider bars)
    geom_text(data = filter(tissue_sig, log10_p >= tissue_threshold),
              aes(label = n_label), 
              hjust = 1.1, size = 3, color = "white", fontface = "bold") +
    # Text outside bar (for narrow bars)
    geom_text(data = filter(tissue_sig, log10_p < tissue_threshold),
              aes(label = n_label), 
              hjust = -0.1, size = 3, color = "black", fontface = "bold") +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05))
    ) +
    coord_flip() +
    labs(
      x = NULL,
      y = expression(-Log[10](italic(P)~value))
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 9, color = "black",
                                 face = tissue_sig$font_face[match(levels(tissue_sig$tissue_order), 
                                                                   tissue_sig$tissue_label)]),
      axis.text.x = element_text(size = 9, color = "black"),
      axis.title.x = element_text(size = 10),
      panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor.x = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      plot.margin = margin(5, 1, 5, 5)
    )
  
  # Right panel: Log2(OR)
  p_tissue_right <- ggplot(tissue_sig, aes(x = tissue_order, y = log2_OR)) +
    geom_col(fill = "steelblue", color = "black", width = 0.85, linewidth = 0.2) +
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
  
  # Combine tissue panels
  combined_tissue <- p_tissue_left + plot_spacer() + p_tissue_right + 
    plot_layout(widths = c(1.3, 0.05, 1)) +
    plot_annotation(
      title = "Tissue Signature Enrichment (TiSigBLUEPRINT)"
    )
  
  print(combined_tissue)
  
  # Save plot
  ggsave("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/Horvath_tissue_signature_enrichment_TiSigBLUEPRINT.png", 
         combined_tissue, width = 10, height = 6, dpi = 300)
  
} else {
  cat("No significant tissue enrichment (TiSigBLUEPRINT) found at p < 0.05 with overlap > 1\n")
  cat("Top tissues by p-value:\n")
  print(tissue_annotated %>% 
          filter(overlap > 0) %>% 
          arrange(p.value) %>% 
          select(tissue_name, overlap, estimate, p.value, FDR) %>%
          head(10))
}

# Save tissue enrichment results
write.csv(tissue_annotated,
          "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/Horvath_tissue_signature_enrichment_TiSigBLUEPRINT_results.csv",
          row.names = FALSE)

# Probe enrichment for CpG context ####
cgi_res <- testEnrichment(query_cpgs, platform="MSA", databases = "KYCG.MSA.CGI.20220904")

# Load CGI database to get context names
load("~/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/MSA_KYCG_Knowledgebases/KYCG.MSA.CGI.20220904.rda")

# Extract CGI context names
cgi_names <- sapply(KYCG.MSA.CGI.20220904, function(x) {
  gsub("CGI;", "", attr(x, "dbname"))
})

# Create mapping dataframe
cgi_mapping <- data.frame(
  cgi_index = as.character(1:length(cgi_names)),
  context_name = cgi_names,
  stringsAsFactors = FALSE
)

# Run enrichment
cgi_res <- testEnrichment(query_cpgs, platform="MSA", databases = "KYCG.MSA.CGI.20220904")

# Add CGI index from rownames
cgi_res$cgi_index <- rownames(cgi_res)

# Prepare CGI enrichment data with proper labels
cgi_annotated <- cgi_res %>%
  as.data.frame() %>%
  left_join(cgi_mapping, by = "cgi_index") %>%
  mutate(
    context = case_when(
      context_name == "Island" ~ "Island",
      context_name == "Shore" ~ "Shore",
      context_name == "Shelf" ~ "Shelf",
      context_name == "OpenSea" ~ "Open Sea",
      TRUE ~ context_name
    ),
    log2_OR = estimate,  # Already in log2 scale
    log10_p = -log10(p.value)
  )

# Create schematic layout with more intuitive colors
schematic_layout <- data.frame(
  context = factor(c("Open Sea", "Shelf", "Shore", "Island"),
                   levels = c("Open Sea", "Shelf", "Shore", "Island")),
  x_start = c(0, 2.5, 5, 7.5),
  x_end = c(2.5, 5, 7.5, 10),
  y_bottom = 0,
  y_top = 1,
  # Colors: Open Sea (dark blue), Shelf (lighter blue), Shore (light blue/cyan), Island (green/yellow)
  region_color = c("#2166ac", "#4393c3", "#92c5de", "#fee090")
)

# Merge with enrichment data
plot_data <- schematic_layout %>%
  left_join(cgi_annotated, by = "context") %>%
  mutate(
    x_mid = (x_start + x_end) / 2,
    significant = !is.na(FDR) & FDR < 0.05,
    label_face = ifelse(significant, "bold", "plain")  # Add font face based on FDR
  )

# Create the plot
cpg_context_plot <- ggplot(plot_data) +
  # Draw rectangles with colors
  geom_rect(aes(xmin = x_start, xmax = x_end, 
                ymin = y_bottom, ymax = y_top, 
                fill = region_color),
            color = "black", linewidth = 1.2) +
  scale_fill_identity() +
  
  # Add significance borders for FDR < 0.05
  geom_rect(data = filter(plot_data, significant),
            aes(xmin = x_start, xmax = x_end, 
                ymin = y_bottom, ymax = y_top),
            color = "red", fill = NA, linewidth = 2.5) +
  
  # Add context labels at top - with conditional bold based on FDR
  geom_text(aes(x = x_mid, y = 1.15, label = context, fontface = label_face), 
            size = 5) +
  
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

print(cpg_context_plot)

ggsave("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/Horvath_cpg_island_context_enrichment.png", cpg_context_plot, width = 10, height = 6, dpi = 300)

write.csv(cgi_annotated, 
          "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/Horvath_cgi_enrichment_results.csv",
          row.names = FALSE)

# Probe genomic proximity of CpGs ####
proximity_res <- testProbeProximity(query_cpgs, platform="MSA")
print(proximity_res)

# Handle inconsistent column naming (P.val vs p.val - changes based on significance!)
p_value <- if ("P.val" %in% names(proximity_res$Stats)) {
  proximity_res$Stats$P.val
} else if ("p.val" %in% names(proximity_res$Stats)) {
  proximity_res$Stats$p.val
} else {
  NA
}

# Check if clustering is significant
# When significant: proximity_res$Clusters is a data frame
# When not significant: proximity_res$Clusters is NA
has_sig_cluster <- !is.na(p_value) && 
  p_value < 0.05 &&
  is.data.frame(proximity_res$Clusters)

# Step 1: Get annotations for CpGs
anno <- sesameData_getManifestGRanges("MSA")
query_anno <- anno[names(anno) %in% query_cpgs]

# Convert to data frame
cpg_positions <- as.data.frame(query_anno) %>%
  mutate(CpG = names(query_anno)) %>%
  select(CpG, seqnames, start)

# Automatically identify cluster members if significant clustering found
if (has_sig_cluster) {
  
  cat("Significant clustering detected!\n")
  
  # Extract cluster information from proximity_res
  cluster_info <- proximity_res$Clusters
  
  cat("\nCluster regions identified:\n")
  print(cluster_info)
  
  # Extract CpGs that fall within cluster regions
  cluster_cpgs <- cpg_positions %>%
    filter(
      # Match chromosome and position range for each cluster
      purrr::pmap_lgl(list(seqnames, start), function(chr, pos) {
        any(
          cluster_info$seqnames == chr &
            pos >= cluster_info$start &
            pos <= cluster_info$end
        )
      })
    )
  
  cat("\nCpGs in the cluster(s):\n")
  print(cluster_cpgs)
  
  # Add cluster membership
  cpg_positions <- cpg_positions %>%
    mutate(in_cluster = CpG %in% cluster_cpgs$CpG)
  
} else {
  # No significant clustering
  cpg_positions <- cpg_positions %>%
    mutate(in_cluster = FALSE)
  
  cluster_cpgs <- data.frame()  # Empty dataframe
  cat("No significant clustering detected (P =", p_value, ")\n")
}

# Step 2: Plot CpG distribution across chromosomes
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

# Create base plot with expanded x-axis for annotation space
cpg_genomic_proximity <- ggplot(cpg_summary, aes(x = seqnames, y = n_cpgs, fill = has_cluster)) +
  geom_col() +
  scale_fill_manual(values = c("FALSE" = "gray60", "TRUE" = "red"),
                    name = "Contains cluster") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(expand = expansion(add = c(3, 3))) +  # Add padding on both sides
  labs(x = "Chromosome", y = "Number of CpGs",
       title = "Distribution of Hypermethylated CpGs Across Chromosomes",
       subtitle = paste0("Spatial clustering P = ", format(p_value, digits = 3))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Only add cluster annotations if significant clustering found
if (has_sig_cluster) {
  
  # Get the chromosome with the cluster for annotation positioning
  cluster_chr <- as.character(cluster_cpgs$seqnames[1])
  chr_num <- gsub("chr", "", cluster_chr)
  
  # Determine if we should annotate to the left or right
  # Right side for chr1-12, left side for chr13-Y
  if (chr_num %in% c(1:12)) {
    annotate_right <- TRUE
    x_offset <- 2  # Offset to the right
  } else {
    annotate_right <- FALSE
    x_offset <- -2  # Offset to the left
  }
  
  # Get numeric chromosome position for positioning
  chr_x <- which(levels(cpg_summary$seqnames) == cluster_chr)
  
  # Prepare cluster annotation data
  cluster_annotation <- cluster_cpgs %>%
    mutate(
      seqnames = factor(seqnames, levels = all_chroms),
      label = CpG
    ) %>%
    left_join(cpg_summary %>% select(seqnames, n_cpgs), by = "seqnames") %>%
    arrange(start) %>%
    mutate(
      label_y = seq(max(cpg_summary$n_cpgs) * 1.1, 
                    max(cpg_summary$n_cpgs) * 1.5, 
                    length.out = dplyr::n())
    )
  
  # Add cluster annotations to plot
  cpg_genomic_proximity <- cpg_genomic_proximity +
    # Add line segments from bar to labels
    geom_segment(data = cluster_annotation,
                 aes(x = seqnames, xend = seqnames, 
                     y = n_cpgs, yend = label_y - 0.5),
                 inherit.aes = FALSE,
                 color = "red", linewidth = 0.5, linetype = "dashed") +
    # Add text labels for CpGs
    geom_text(data = cluster_annotation,
              aes(x = seqnames, y = label_y, label = label),
              inherit.aes = FALSE,
              size = 3, fontface = "bold", color = "red")
  
  # Add distance annotation if exactly 2 CpGs in cluster
  if (nrow(cluster_annotation) == 2) {
    distance_bp <- abs(diff(cluster_annotation$start))
    distance_y <- mean(cluster_annotation$label_y)
    distance_label <- paste0(distance_bp, " bp")
    
    if (annotate_right) {
      # Annotate to the RIGHT (chr1-12)
      cpg_genomic_proximity <- cpg_genomic_proximity +
        # Add distance label to the right
        annotate("text", 
                 x = chr_x + x_offset, 
                 y = distance_y,
                 label = distance_label,
                 size = 3.5, fontface = "italic", color = "darkred", hjust = 0) +
        # Add arrow to UPPER CpG from the right
        annotate("segment",
                 x = chr_x + x_offset - 0.1,
                 xend = chr_x,
                 y = distance_y,
                 yend = cluster_annotation$label_y[2] - 0.3,
                 color = "darkred", linewidth = 0.6, 
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
        # Add arrow to LOWER CpG from the right
        annotate("segment",
                 x = chr_x + x_offset - 0.1,
                 xend = chr_x,
                 y = distance_y,
                 yend = cluster_annotation$label_y[1] + 0.3,
                 color = "darkred", linewidth = 0.6,
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed"))
    } else {
      # Annotate to the LEFT (chr13-Y)
      cpg_genomic_proximity <- cpg_genomic_proximity +
        # Add distance label to the left
        annotate("text", 
                 x = chr_x + x_offset, 
                 y = distance_y,
                 label = distance_label,
                 size = 3.5, fontface = "italic", color = "darkred", hjust = 1) +
        # Add arrow to UPPER CpG from the left
        annotate("segment",
                 x = chr_x + x_offset + 0.1,
                 xend = chr_x,
                 y = distance_y,
                 yend = cluster_annotation$label_y[2] - 0.3,
                 color = "darkred", linewidth = 0.6, 
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
        # Add arrow to LOWER CpG from the left
        annotate("segment",
                 x = chr_x + x_offset + 0.1,
                 xend = chr_x,
                 y = distance_y,
                 yend = cluster_annotation$label_y[1] + 0.3,
                 color = "darkred", linewidth = 0.6,
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed"))
    }
  }
}

# Step 3: Get genes near the cluster CpGs
if (has_sig_cluster) {
  gene_dbs <- buildGeneDBs(cluster_cpgs$CpG, platform = "MSA")
  gene_res <- testEnrichment(cluster_cpgs$CpG, gene_dbs, platform = "MSA")
  
  cluster_genes <- gene_res %>%
    as.data.frame() %>%
    filter(overlap > 0) %>%
    arrange(desc(overlap))
  
  cat("\nGenes associated with cluster CpGs:\n")
  print(cluster_genes %>% select(gene_name, overlap, estimate, p.value))
}

print(cpg_genomic_proximity)

ggsave("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/DunedinPACE_cpg_genomic_proximity.png", 
       cpg_genomic_proximity, 
       width = 10, 
       height = 6, 
       dpi = 300)

# Save the test statistics
proximity_stats <- data.frame(
  test = "Genomic Proximity",
  p_value = p_value,
  n_cpgs = length(query_cpgs),
  significant_clustering = has_sig_cluster
)

write.csv(proximity_stats, 
          "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/DunedinPACE_proximity_test_results.csv",
          row.names = FALSE)

# Save chromosome distribution
write.csv(cpg_summary, 
          "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/DunedinPACE_cpg_chromosome_distribution.csv",
          row.names = FALSE)

# If significant clustering found, save cluster details
if (has_sig_cluster) {
  write.csv(cluster_cpgs, 
            "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/DunedinPACE_cpg_cluster_details.csv",
            row.names = FALSE)
  
  # Save cluster-associated genes if they exist
  if (exists("cluster_genes")) {
    write.csv(cluster_genes, 
              "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/DunedinPACE_cluster_associated_genes.csv",
              row.names = FALSE)
  }
}

# Probe enrichment for genomic region ####
# Metagene enrichment analysis
metagene_res <- testEnrichment(
  query_cpgs,
  platform = "MSA",
  databases = "KYCG.MSA.MetagenePC.20220911"
)

# Load the metagene database to get labels
load("~/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/MSA_KYCG_Knowledgebases/KYCG.MSA.MetagenePC.20220911.rda")

# Extract metagene region names and parse them keeping prefix for uniqueness
metagene_names <- sapply(KYCG.MSA.MetagenePC.20220911, function(x) {
  full_name <- gsub("MetagenePC;", "", attr(x, "dbname"))
  full_name  # Keep the full label including prefix for now
})

# Create mapping dataframe with better formatted labels
metagene_mapping <- data.frame(
  metagene_index = as.character(1:length(metagene_names)),
  region_code = metagene_names,
  stringsAsFactors = FALSE
) %>%
  mutate(
    # Parse prefix and suffix
    prefix = as.numeric(sub(";.*", "", region_code)),
    suffix = sub(".*;", "", region_code),
    # Create clean labels - keep prefix only for non-percentage regions to distinguish them
    region_label = ifelse(grepl("%", suffix), 
                          suffix,  # Just use percentage for gene body
                          paste0("[", prefix, "] ", suffix))  # Add prefix in brackets for distance regions
  )

# Add metagene index from rownames
metagene_res$metagene_index <- rownames(metagene_res)

# Add labels to metagene_df
metagene_df <- metagene_res %>%
  as.data.frame() %>%
  left_join(metagene_mapping, by = "metagene_index") %>%
  arrange(prefix) %>%  # Order by prefix number (which goes from -10 to 19)
  mutate(
    region_order = factor(region_label, levels = region_label),  # Preserve sorted order
    log2_OR = estimate  # Already in log2 scale
  )

# Display results
print(metagene_df %>% select(region_label, estimate, overlap, p.value, FDR))

# Check if any regions have FDR < 0.05
has_fdr_sig <- any(metagene_df$FDR < 0.05, na.rm = TRUE)

# Determine which y-axis to use
if (has_fdr_sig) {
  # Use FDR if there are significant results
  metagene_df <- metagene_df %>%
    mutate(y_value = -log10(FDR))
  
  y_label <- expression(-Log[10](FDR))
  hline_value <- -log10(0.05)
  subtitle_text <- "Location relative to gene structure"
  
} else {
  # Use p-value if no FDR significant results
  metagene_df <- metagene_df %>%
    mutate(y_value = -log10(p.value))
  
  y_label <- expression(-Log[10](italic(P)~value))
  hline_value <- -log10(0.05)
  subtitle_text <- "Location relative to gene structure (No significant enrichment after FDR correction)"
}

# Visualize gene regions with conditional y-axis
p_metagene <- ggplot(metagene_df, aes(x = region_order, y = y_value)) +
  geom_col(fill = "gray60", color = "black", width = 0.7) +
  geom_text(aes(label = overlap), vjust = -0.5, size = 3) +
  geom_hline(yintercept = hline_value, linetype = "dashed", color = "red") +
  labs(
    title = "Gene Region Enrichment",
    subtitle = subtitle_text,
    x = "Gene Region",
    y = y_label
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

print(p_metagene)

ggsave("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/Horvath_enrichment_for_genomic_region.png", p_metagene, width = 10, height = 6, dpi = 300)

write.csv(metagene_df, 
          "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/Horvath_metagene_enrichment_results.csv",
          row.names = FALSE)
