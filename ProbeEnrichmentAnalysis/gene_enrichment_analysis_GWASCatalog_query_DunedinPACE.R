# This script runs a gene enrichment analysis, using KnowYourCG, for DunedinPACE CpGs that are differently methylated by genetic ancestry.
# Enriched genes are then queried against the GWAS Catalog to identify links to neurodegenerative or cognitive aging traits.
# Author: Rory Boyle rorytboyle@gmail.com
# Date: 14/11/2025
# Updated: 2nd December 2025 to use development version of KnowYourCG and use MSA platform
# Updated: 9th December 2025 to update cpg list and to fix sig_genes_plot (to plot top 25 enriched genes instead of all)

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
library(patchwork)
library(gwasrapidd)
library(gprofiler2)

# Extract differentially methylated CpGs ####
# Cache databases (only needed once)
sesameDataCache()

# Load your differential methylation results
query_cpgs <- readRDS("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251208_hypermethylated_DunedinPACE_CpGs_African_Ancestry.rds")

cat(sprintf("Analyzing %d significant CpGs\n", length(query_cpgs)))

# Gene enrichment analysis ####

# Map CpGs to genes and test enrichment
gene_dbs <- buildGeneDBs(query_cpgs, platform = "MSA")
gene_res <- testEnrichment(query_cpgs, gene_dbs, platform = "MSA")

write.csv(gene_res, '/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/DunedinPACE_gene_enrichment_results.csv', row.names = FALSE)

# Get significant genes
sig_genes <- gene_res %>%
  as.data.frame() %>%
  filter(FDR < 0.05) %>%
  arrange(FDR)

cat(sprintf("Found %d significantly enriched genes (FDR < 0.05)\n", nrow(sig_genes)))

## Make plot for all FDR sig genes ####
sig_genes_plot <- sig_genes %>%
  # Filter for top 25
  slice_max(order_by = -log10.p.value, n = 25) %>%
    mutate(
    gene_label = gene_name, 
    n_label = paste0("N = ", overlap),
    log2_OR = log2(estimate),
    log10_p = -log10(p.value),
    gene_order = factor(gene_label, levels = rev(gene_label))
  )

cat(sprintf("All FDR significant genes: %d\n", nrow(sig_genes_plot)))

# Left panel with smaller text
p_gene_left <- ggplot(sig_genes_plot, aes(x = gene_order, y = log10_p)) +
  geom_col(fill = "grey40", color = "black", width = 0.85, linewidth = 0.2) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_text(aes(label = n_label), 
            hjust = 1.1,
            size = 2,  # Reduced from 3
            color = "white",
            fontface = "bold") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)), 
    breaks = c(0, 2, 4, 6, 8), 
    labels = c("1", expression(10^-2), expression(10^-4), expression(10^-6), expression(10^-8)) 
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = expression(-Log[10](italic(P)~value))
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 6, color = "black", face = "bold"),  # Reduced from 9
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 10),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    plot.margin = margin(5, 1, 5, 5)
  ) 

# Right panel with smaller text
p_gene_right <- ggplot(sig_genes_plot, aes(x = gene_order, y = log2_OR)) + 
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

# Combine and save with increased height
combined_genes <- p_gene_left + plot_spacer() + p_gene_right + 
  plot_layout(widths = c(1.3, 0.05, 1)) +
  plot_annotation(
    title = "Probe enrichment for genes (top 25 by -Log10 p-value)"
  )

ggsave("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/DunedinPACE_gene_enrichment_analysis_overlap_all.png", 
       combined_genes, 
       width = 6, 
       height = 14,  # Tall figure
       dpi = 300)

combined_genes

## Make plot for only genes with overlap > 1 ####
sig_genes_plot_overlap_greaterThan1 <- sig_genes %>%
  # Filter 1): Keep only genes hit by more than one significant CpG
  filter(overlap > 1) %>%
  # Prepare plotting variables
  mutate(
    gene_label = gene_name, 
    n_label = paste0("N = ", overlap),
    log2_OR = log2(estimate),
    log10_p = -log10(p.value),
    gene_order = factor(gene_label, levels = rev(gene_label))
  )

cat(sprintf("Genes with >1 CpG overlap: %d\n", nrow(sig_genes_plot_overlap_greaterThan1)))

# Left panel: P-values with N labels embedded in bars
p_gene_left_overlap_greaterThan1 <- ggplot(sig_genes_plot_overlap_greaterThan1, aes(x = gene_order, y = log10_p)) +
  geom_col(fill = "grey40", color = "black", width = 0.85, linewidth = 0.2) +
  # add red dashed line to indicate nominal singificance
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_text(aes(label = n_label), 
            hjust = 1.1,
            size = 3,
            color = "white",
            fontface = "bold") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)), 
    breaks = c(0, 2, 4, 6, 8), 
    labels = c("1", expression(10^-2), expression(10^-4), expression(10^-6), expression(10^-8)) 
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = expression(-Log[10](italic(P)~value))
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 9, color = "black", face = "bold"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 10),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    plot.margin = margin(5, 1, 5, 5)
  ) 

# Right panel: Log2(OR)
p_gene_right_overlap_greaterThan1 <- ggplot(sig_genes_plot_overlap_greaterThan1, aes(x = gene_order, y = log2_OR)) + 
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
combined_genes_overlap_greaterThan1 <- p_gene_left_overlap_greaterThan1 + plot_spacer() + p_gene_right_overlap_greaterThan1 + 
  plot_layout(widths = c(1.3, 0.05, 1)) +
  plot_annotation(
    title = "Probe enrichment for genes"
  )

print(combined_genes_overlap_greaterThan1) # is filtering to overlap = 2 ok?

# Save
ggsave("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/DunedinPACE_gene_enrichment_analysis_overlap_greaterThan1.png", combined_genes_overlap_greaterThan1, width = 6, height = 4, dpi = 300)

# Query GWAS catalog for enriched genes ####
my_genes <- sig_genes %>%
  # Filter 1): Keep only genes hit by more than one significant CpG
  filter(overlap > 1)

my_genes <- my_genes$gene_name

# Comprehensive list for neurodegenerative/cognitive aging
neurodegenerative_traits <- c(
  # Core cognitive terms
  "cognitive",
  "cognition",
  "memory",
  "executive function",
  
  # Dementia types
  "alzheimer",
  "dementia",
  "lewy body",
  "frontotemporal",
  "vascular dementia",
  
  # Parkinson's related
  "parkinson",
  
  # Other
  # "PART", # this pulls 1 EFO_0004736   aspartate aminotransferase measurement 
  "primary age-related tauopathy",
  "als",
  
  # Related conditions
  "mild cognitive",
  "neurodegener",
  "brain aging",
  "cognitive aging",
  "MCI",
  "brain age",
  
  # Biomarkers
  "amyloid",
  "amyloid-beta", 
  "tau",
  "tdp-43",
  "alpha-synuclein",
  "neurofilament light chain",
  "GFAP",
  "neuroimaging",
  
  # Neuropathology
  "neurofibrillary tangles",
  "neuritic plaques",
  
  # Related phenotypes
  "hippocampal volume",
  "cortical thickness",
  "white matter",
  "brain atrophy"
)

# Create pattern for neurodegenerative traits
trait_pattern <- paste(neurodegenerative_traits, collapse = "|")

# Initialize results
all_results <- data.frame()

# For each gene, get all associations, then filter for neurodegenerative traits
for (gene in my_genes) {
  cat("Querying", gene, "...\n")
  
  # Get variants for this gene
  variants <- get_variants(gene_name = gene)
  
  if (gwasrapidd::n(variants) > 0) {
    # Get associations for these variants
    assocs <- get_associations(variant_id = variants@variants$variant_id)
    
    if (gwasrapidd::n(assocs) > 0) {
      # Get studies for these associations
      studies <- get_studies(association_id = assocs@associations$association_id)
      
      # Get traits for these associations
      traits <- get_traits(association_id = assocs@associations$association_id)
      
      # Filter studies for neurodegenerative reported traits
      neuro_reported <- studies@studies %>%
        filter(grepl(trait_pattern, reported_trait, ignore.case = TRUE)) # where does reported trait come from here? within studies@studies?
      
      # Filter traits for neurodegenerative EFO traits
      neuro_efo <- traits@traits %>%
        filter(grepl(trait_pattern, trait, ignore.case = TRUE)) # where does EFO trait come from here? within traits@traits?
      
      # Build table for reported traits
      if (nrow(neuro_reported) > 0) {
        for (j in 1:nrow(neuro_reported)) {
          study_id <- neuro_reported$study_id[j]
          reported_trait <- neuro_reported$reported_trait[j]
          
          # Get associations for this study
          study_assocs <- get_associations(study_id = study_id)
          
          if (gwasrapidd::n(study_assocs) > 0) {
            # Get all association IDs from this study
            study_assoc_ids <- study_assocs@associations$association_id
            
            # Get associations for our gene
            gene_variants <- get_variants(gene_name = gene)
            if (gwasrapidd::n(gene_variants) > 0) {
              gene_all_assocs <- get_associations(variant_id = gene_variants@variants$variant_id)
              gene_assoc_ids <- gene_all_assocs@associations$association_id
              
              # Find overlap between study associations and gene associations
              matching_assoc_ids <- intersect(study_assoc_ids, gene_assoc_ids)
              
              if (length(matching_assoc_ids) > 0) {
                for (assoc_id in matching_assoc_ids) {
                  assoc_details <- study_assocs@associations %>%
                    filter(association_id == assoc_id) %>%
                    slice(1)
                  
                  if (nrow(assoc_details) > 0) {
                    all_results <- rbind(all_results, data.frame(
                      gene = gene,
                      trait = reported_trait,
                      trait_type = "reported",
                      association_id = assoc_id,
                      study_id = study_id,
                      pvalue = assoc_details$pvalue[1],
                      or_per_copy = ifelse(!is.null(assoc_details$or_per_copy_number[1]) && 
                                             !is.na(assoc_details$or_per_copy_number[1]), 
                                           assoc_details$or_per_copy_number[1], NA),
                      beta = ifelse(!is.null(assoc_details$beta_number[1]) && 
                                      !is.na(assoc_details$beta_number[1]), 
                                    assoc_details$beta_number[1], NA),
                      study_url = paste0("https://www.ebi.ac.uk/gwas/studies/", study_id),
                      stringsAsFactors = FALSE
                    ))
                  }
                }
              }
            }
          }
        }
      }
      
      # Build table for EFO traits
      if (nrow(neuro_efo) > 0) {
        for (k in 1:nrow(neuro_efo)) {
          efo_id <- neuro_efo$efo_id[k]
          efo_trait <- neuro_efo$trait[k]
          
          # Get associations for this EFO trait
          efo_assocs <- get_associations(efo_id = efo_id)
          
          if (gwasrapidd::n(efo_assocs) > 0) {
            # Get all association IDs from this EFO trait
            efo_assoc_ids <- efo_assocs@associations$association_id
            
            # Check which of these associations involve our gene
            # Use the original assocs or query the gene directly
            gene_variants <- get_variants(gene_name = gene)
            if (gwasrapidd::n(gene_variants) > 0) {
              gene_all_assocs <- get_associations(variant_id = gene_variants@variants$variant_id)
              gene_assoc_ids <- gene_all_assocs@associations$association_id
              
              # Find overlap between EFO associations and gene associations
              matching_assoc_ids <- intersect(efo_assoc_ids, gene_assoc_ids)
              
              if (length(matching_assoc_ids) > 0) {
                for (assoc_id in matching_assoc_ids) {
                  assoc_details <- efo_assocs@associations %>%
                    filter(association_id == assoc_id) %>%
                    slice(1)
                  
                  if (nrow(assoc_details) > 0) {
                    assoc_study <- get_studies(association_id = assoc_id)
                    
                    if (gwasrapidd::n(assoc_study) > 0) {
                      all_results <- rbind(all_results, data.frame(
                        gene = gene,
                        trait = efo_trait,
                        trait_type = "EFO",
                        association_id = assoc_id,
                        study_id = assoc_study@studies$study_id[1],
                        pvalue = assoc_details$pvalue[1],
                        or_per_copy = ifelse(!is.null(assoc_details$or_per_copy_number[1]) && 
                                               !is.na(assoc_details$or_per_copy_number[1]), 
                                             assoc_details$or_per_copy_number[1], NA),
                        beta = ifelse(!is.null(assoc_details$beta_number[1]) && 
                                        !is.na(assoc_details$beta_number[1]), 
                                      assoc_details$beta_number[1], NA),
                        study_url = paste0("https://www.ebi.ac.uk/gwas/studies/", 
                                           assoc_study@studies$study_id[1]),
                        stringsAsFactors = FALSE
                      ))
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

# Add specific variants to each of the gene-trait associations
all_results_with_variants <- data.frame()
for (i in 1:nrow(all_results)) {
  assoc_id <- all_results$association_id[i]
  
  # Get association details
  assoc_details <- get_associations(association_id = assoc_id)
  
  # Initialize variant info
  variant_ids <- NA
  risk_allele_info <- NA
  risk_freq <- NA
  
  if (gwasrapidd::n(assoc_details) > 0) {
    # Get risk alleles
    risk_alleles <- assoc_details@risk_alleles
    
    if (nrow(risk_alleles) > 0) {
      # Get variant IDs
      variant_ids <- paste(risk_alleles$variant_id, collapse = "; ")
      
      # Get risk alleles
      risk_allele_info <- paste(risk_alleles$risk_allele, collapse = "; ")
      
      # Get risk frequencies if available
      if ("risk_frequency" %in% colnames(risk_alleles)) {
        risk_freq <- paste(risk_alleles$risk_frequency, collapse = "; ")
      }
    }
  }
  
  # Add to results
  all_results_with_variants <- rbind(all_results_with_variants, data.frame(
    all_results[i, ],
    variant_id = variant_ids,
    risk_allele = risk_allele_info,
    risk_frequency = risk_freq,
    stringsAsFactors = FALSE
  ))
  
  if (i %% 10 == 0) cat(".")
}

# Clean final_results 
final_results <- all_results_with_variants %>%
  # remove non-relevant traits
  filter(!trait %in% c("LGALS3BP protein levels", # sneaks in because of search for ALS
                       "Cognitive empathy")) %>% # sneaks in because of search for cognitive (cognitive empathy more of a psychological construct)
  arrange(association_id, desc(trait_type)) %>%
  distinct(association_id, .keep_all = TRUE) %>%
  arrange(gene, pvalue)


cat("Found", nrow(final_results), "variant-trait associations involving these genes:", unique(final_results$gene), "\n\n")

# Create a version with a single gene-trait association
final_results_unique <- final_results %>%
  arrange(gene, trait, pvalue) %>%
  distinct(gene, trait, .keep_all = TRUE) %>%
  arrange(gene, pvalue)

cat("Found", nrow(final_results_unique), "unique variant-trait associations involving these genes:", unique(final_results$gene), "\n\n")

# Plot gene-trait associations
# Prepare data with aggregation
viz_data <- final_results %>%
  mutate(
    trait_wrapped = str_wrap(trait, width = 40),
    variant_label = variant_id
  ) %>%
  arrange(gene, variant_id, trait_wrapped) %>%
  # Group by position and aggregate
  group_by(gene, variant_label, trait_wrapped) %>%
  summarize(
    n_associations = dplyr::n(),
    best_pvalue = min(pvalue, na.rm = TRUE),
    neg_log10_p = -log10(best_pvalue),
    .groups = "drop"
  ) %>%
  mutate(trait_wrapped = factor(trait_wrapped, levels = unique(trait_wrapped)))

p <- ggplot(viz_data, aes(x = variant_label, 
                          y = trait_wrapped,
                          size = n_associations, 
                          color = neg_log10_p)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(2, 10), name = "N Associations", 
                        breaks = function(x) unique(floor(pretty(seq(min(x), max(x)))))) +
  scale_color_gradient(low = "lightblue", high = "darkred",
                       name = expression(-log[10](italic(P)))) +
  facet_grid(. ~ gene, scales = "free_x", space = "free_x") +
  labs(x = "Variant", y = "Trait", 
       title = "Gene-Trait Associations from GWAS Catalog") +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 8, lineheight = 0.9),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray90", color = "gray60"),
    legend.position = "right",
    panel.spacing = unit(0.5, "lines")
  )


print(p)
# Save
write.csv(final_results, "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/DunedinPACE_gene_enrichment_overlap_greaterThan1_neurodegenerative_associations.csv", row.names = FALSE)
ggsave("/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/DunedinPACE_gene_enrichment_overlap_greaterThan1_neurodegenerative_associations_plot.png",
       p, 
       width = 10.3, 
       height = 5.47, 
       units = "in",
       dpi = 300)