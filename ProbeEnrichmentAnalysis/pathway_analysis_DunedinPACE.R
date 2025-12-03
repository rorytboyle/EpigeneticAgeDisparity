# This script runs a pathway analysis, using gprofiler2, for genes that are enriched for CpGs that are differently methylated by genetic ancestry.
# Author: Rory Boyle rorytboyle@gmail.com
# Date: 17/11/2025

library(gprofiler2)

# Load gene enrichment analysis results ####
# from gene_enrichment_analysis_with_GWASCatalog.R
gene_res <- read_csv('/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/DunedinPACE_gene_enrichment_results.csv')

# Get significant genes (FDR < 0.05), hit by more than one significant CpG
my_genes <- gene_res %>%
  as.data.frame() %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  # Filter 1): Keep only genes hit by more than one significant CpG
  filter(overlap > 1)

my_genes <- my_genes$gene_name

# Define file paths for gene annotations if already downloaded
gene_annot_file <- "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/data/msa_genes_annot_hg38.rds"
background_genes_file <- "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/data/epic_background_genes_hg38.rds"

# Get gene annotation and background gene list ####
# Download and process gene annotation if not already cached
if (file.exists(background_genes_file)) {
  cat("Loading cached background genes from:", background_genes_file, "\n")
  all_epic_genes_hg38 <- readRDS(background_genes_file)
  cat(sprintf("Loaded %d unique genes\n", length(all_epic_genes_hg38)))
} else {
  cat("Background genes not found, downloading MSA gene annotation...\n")
  
  # Download MSA gene annotation (hg38, GENCODE v41)
  gene_annot_url <- "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/MSA/MSA.hg38.manifest.gencode.v41.tsv.gz"
  
  temp_file <- tempfile(fileext = ".tsv.gz")
  download.file(gene_annot_url, temp_file, mode = "wb", quiet = TRUE)
  
  # Read annotation
  msa_genes_annot <- read.delim(temp_file, header = TRUE)
  
  cat("Columns in gene annotation file:\n")
  print(colnames(msa_genes_annot))
  
  cat("\nFirst few rows:\n")
  print(head(msa_genes_annot))
  
  # Extract gene names
  all_epic_genes_hg38 <- msa_genes_annot$geneNames
  
  # Split if multiple genes per probe (likely semicolon-separated)
  all_epic_genes_hg38 <- unique(unlist(strsplit(as.character(all_epic_genes_hg38), ";")))
  all_epic_genes_hg38 <- trimws(all_epic_genes_hg38)  # Remove whitespace
  all_epic_genes_hg38 <- all_epic_genes_hg38[!is.na(all_epic_genes_hg38) & 
                                               all_epic_genes_hg38 != "" &
                                               all_epic_genes_hg38 != "NA"]
  
  cat(sprintf("\nBackground: %d unique genes on MSA/EPIC array (hg38, GENCODE v41)\n", 
              length(all_epic_genes_hg38)))
  
  # Save for future use
  saveRDS(msa_genes_annot, gene_annot_file)
  saveRDS(all_epic_genes_hg38, background_genes_file)
  cat("\nSaved background gene list to:", background_genes_file, "\n")
  
  # Clean up
  unlink(temp_file)
}

# Run pathway enrichment ####
cat(sprintf("\nRunning pathway analysis on %d genes...\n", length(my_genes)))

pathway_res <- gost(
  query = my_genes,
  organism = "hsapiens",
  sources = c("GO:BP", "KEGG", "REAC", "WP"),
  custom_bg = all_epic_genes_hg38,
  correction_method = "fdr",
  significant = FALSE,  # Only return significant results
  user_threshold = 0.05,
  evcodes = TRUE
)

# Process results
if (!is.null(pathway_res$result)) {
  
  sig_pathways <- pathway_res$result %>%
    filter(p_value < 0.05) %>%
    # filter(term_size >= 10 & term_size <= 500) %>%  # Remove overly broad/specific terms
    mutate(
      # gene_ratio = intersection_size / length(my_genes),
      # bg_ratio = term_size / length(all_epic_genes_hg38),
      # enrichment_fold = gene_ratio / bg_ratio,
      # neg_log10_p = -log10(p_value)
    ) %>%
    arrange(p_value)
  
  cat(sprintf("\nFound %d significant pathways (p < 0.05, term size 10-500)\n", 
              nrow(sig_pathways)))
  
  if (nrow(sig_pathways) > 0) {
    # Display top pathways
    cat("\nTop 20 enriched pathways:\n")
    print(sig_pathways %>%
            select(source, term_name, p_value, intersection_size, term_size, 
                   ) %>%
            head(20))

    
  } else {
    cat("\nNo pathways passed filtering criteria (p < 0.05, term size 10-500)\n")
  }
  
} else {
  cat("\nNo significant pathways found\n")
}

# Save output ####
output_rds <- "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/DunedinPACE_pathway_enrichment_results_hg38.rds"
saveRDS(sig_pathways, output_rds)

# Save simplified version as CSV (for easy viewing)
output_csv <- "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/DunedinPACE_pathway_enrichment_results_hg38.csv"
sig_pathways_csv <- sig_pathways %>%
  mutate(
    # Convert gene list to semicolon-separated string
    genes_in_pathway = sapply(intersection, function(x) {
      if (is.list(x)) paste(unlist(x), collapse = ";") else paste(x, collapse = ";")
    })
  ) %>%
  select(-intersection, -any_of("parents"))  # Remove remaining list columns

write.csv(sig_pathways_csv, output_csv, row.names = FALSE)

# Visualize pathway analysis results ####
# make an interactive plot (hover over each point)
p_pathway <- gostplot(pathway_res, capped = FALSE, interactive = TRUE)
print(p_pathway)

# Create table of sig pathways ####
# Create a table of significant pathways
# Top 20 pathways, clean formatting
top20_ids <- pathway_res$result %>%
  arrange(p_value) %>%
  head(20) %>%
  pull(term_id)

pathway_res_top20 <- pathway_res
pathway_res_top20$result <- pathway_res$result %>%
  filter(term_id %in% top20_ids) 

# Create clean table
pathway_table <- publish_gosttable(
  pathway_res_top20,
  use_colors = FALSE,  # Black and white for publication
  show_columns = c("source", "term_name", "p_value", 
                   "intersection_size", "term_size"))

# View the table
print(pathway_table)
ggsave(pathway_table, 
       filename = "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/DunedinPACE_top20_enriched_pathways_hg38.png", 
       width = 8, height = 6, dpi = 300)

# Examine specific pathway (Fatty acid biosynthesis as it is role of ECH...) ####
# # Highlight fatty acid synthesis
# # Highlight specific pathway(s) by their term ID
# # Find the term ID for your pathway of interest first
# # pathway_of_interest <- sig_pathways %>%
# #   filter(grepl("fatty acid biosynthesis", term_name, ignore.case = TRUE))
# 
# # View the term_id
# # print(pathway_of_interest$term_id)
# 
# # Visualize pathways
# # p_pathway_fatty_acid <- gostplot(pathway_res, capped = F, interactive = F)
# # print(p_pathway_fatty_acid)
# 
# # # Highlight it on the plot
# # p_highlighted <- publish_gostplot(
# #   p_pathway_fatty_acid,
# #   highlight_terms = pathway_of_interest$term_id,  # Can be a vector of multiple IDs
# #   width = 8,
# #   height = 6
# # )
# 
# # print(p_highlighted)
# 
# # # Save
# # ggsave(p_highlighted, 
# #        filename = "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/DunedinPACE_enriched_pathways_highlight_fatty_acid_hg38.png", 
# #        width = 8, height = 6, dpi = 300)
