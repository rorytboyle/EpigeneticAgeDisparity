# This script prepares a list of MSA CpGs from differential methylation analyses for enrichment analyses using KnowYourCG
# Author: Rory Boyle rorytboyle@gmail.com
# Date: 26/11/2025

# Run functions from mQTLAnalysis/get_windows.R ####
# Take get_msa_annotation() and lookup_multiple_cpgs() functions from mQTLAnalysis/get_windows.R
# download and load MSA hg38 manifest with genomic positions from Zhou lab github
get_msa_annotation <- function() {
  
  # Download MSA manifest with mapping information (hg38)
  msa_url <- "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/MSA/MSA.hg38.manifest.tsv.gz"
  
  cat("Downloading MSA manifest file...\n")
  temp_file <- tempfile(fileext = ".tsv.gz")
  
  # Download the file
  download.file(msa_url, temp_file, mode = "wb", quiet = TRUE)
  
  cat("Loading annotation data...\n")
  # Read the annotation file
  msa_manifest <- read.table(temp_file, header = TRUE, sep = "\t", 
                             stringsAsFactors = FALSE, quote = "")
  
  # Clean up temp file
  unlink(temp_file)
  
  cat(paste("Loaded", nrow(msa_manifest), "probes from MSA manifest\n"))
  
  return(msa_manifest)
}

manifest <- get_msa_annotation()

# lookup multiple CpGs at once (handles _BC21, _TC21, _TC11, etc. suffixes)
lookup_multiple_cpgs <- function(cpg_list, manifest = NULL) {
  
  # Load manifest if not provided
  if (is.null(manifest)) {
    manifest <- get_msa_annotation()
  }
  
  cat(paste("Looking up", length(cpg_list), "CpGs...\n"))
  
  # Extract base CpG names from manifest (remove suffixes like _TC21, _BC11)
  manifest$Base_CpG <- gsub("_.*$", "", manifest$Probe_ID)
  
  # Initialize results
  all_results <- data.frame()
  found_inputs <- c()
  
  for (cpg_input in cpg_list) {
    # First try exact match
    exact_match <- manifest[manifest$Probe_ID == cpg_input, ]
    
    if (nrow(exact_match) > 0) {
      # Found exact match
      result <- data.frame(
        Input_CpG = cpg_input,
        Manifest_ID = exact_match$Probe_ID[1],
        Chromosome = exact_match$CpG_chrm[1],
        Start = exact_match$CpG_beg[1],
        End = exact_match$CpG_end[1],
        Match_Type = "exact",
        stringsAsFactors = FALSE
      )
      all_results <- rbind(all_results, result)
      found_inputs <- c(found_inputs, cpg_input)
      cat(paste("  Found exact match for:", cpg_input, "->", exact_match$Probe_ID[1], "\n"))
      
    } else {
      # Try base name matching
      base_matches <- manifest[manifest$Base_CpG == cpg_input, ]
      
      if (nrow(base_matches) > 0) {
        # Take the first match if multiple exist
        result <- data.frame(
          Input_CpG = cpg_input,
          Manifest_ID = base_matches$Probe_ID[1],
          Chromosome = base_matches$CpG_chrm[1],
          Start = base_matches$CpG_beg[1],
          End = base_matches$CpG_end[1],
          Match_Type = "base_name",
          stringsAsFactors = FALSE
        )
        all_results <- rbind(all_results, result)
        found_inputs <- c(found_inputs, cpg_input)
        cat(paste("  Found base match for:", cpg_input, "->", base_matches$Probe_ID[1], "\n"))
        
        # Show all available matches if multiple exist
        if (nrow(base_matches) > 1) {
          cat(paste("    Note:", nrow(base_matches), "total matches found for", cpg_input, ":", 
                    paste(base_matches$Probe_ID[1:min(5, nrow(base_matches))], collapse=", "), 
                    if(nrow(base_matches) > 5) "..." else "", "\n"))
        }
      } else {
        cat(paste("  No match found for:", cpg_input, "\n"))
      }
    }
  }
  
  if (nrow(all_results) == 0) {
    cat("No CpGs found in MSA manifest\n")
    return(NULL)
  }
  
  # Add standard columns
  all_results$Genome_Build <- "hg38"
  all_results$Array_Type <- "MSA"
  
  # Report summary
  missing_cpgs <- setdiff(cpg_list, found_inputs)
  if (length(missing_cpgs) > 0) {
    cat(paste("CpGs not found:", paste(missing_cpgs, collapse = ", "), "\n"))
  }
  
  cat(paste("Successfully found", nrow(all_results), "out of", length(cpg_list), "CpGs\n"))
  
  return(all_results)
}

# Prep query (get CpG list) ####
# only select CpGs that are significantly hypermethylated in African Ancestry at FDR < 0.05 in unadjusted and adjusted (for cell type proportion) analyses
# These are the blue dots in the Δβ Change of Clock CpGs by Genetic Ancestry after Cell Type Adjustment plot
dmp_cpgs <- read.csv('/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20251126_DunedinPACE_adjusted_vs_unadjusted_comparison.csv') %>%
  filter(Ancestry_unadj == "Higher in AFR" & Ancestry_adj == "Higher in AFR") %>%
  filter(adj.P.Val_unadj < 0.05 & adj.P.Val_adj < 0.05)

# Get MSA manifest IDs ####
# Match CpGs to MSA format (with suffix like _TC21, _BC11, etc.) from manifest
positions <- lookup_multiple_cpgs(dmp_cpgs$CpG, manifest)

# Retain CpG ids for query for CpG IDs only
query_cpgs <- positions$Manifest_ID

# Save
saveRDS(query_cpgs, "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/results/20250201_hypermethylated_DunedinPACE_CpGs_African_Ancestry.rds")

