# This script identifies windows within a specified distance of a genomic position 
# of provided CpGs from the MSA array for mQTL analyses.
# Author: Rory Boyle rorytboyle@gmail.com
# Date: 09/09/2025

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

# create analysis windows around CpGs for mQTL analysis
create_mqtl_windows <- function(cpg_positions, window_distance = 3000) {
  
  if (is.null(cpg_positions) || nrow(cpg_positions) == 0) {
    cat("No CpG positions provided\n")
    return(NULL)
  }
  
  cat(paste("Creating +/-", window_distance, "bp windows around", nrow(cpg_positions), "CpGs...\n"))
  cat(paste("Total window size:", 2 * window_distance, "bp\n")) # Standard in mQTL studies is to select a distance and look at the CpG location +/- that distance 
                                                                # e.g. if window = 3Kbp look at position -3000 to +3000 around the CpG not -1500 to +1500
  
  # Use CpG_beg as the center point (0-based coordinates)
  cpg_center <- cpg_positions$Start
  
  # Create windows: +/- window_distance from center
  windows <- data.frame(
    Input_CpG = cpg_positions$Input_CpG,
    Manifest_ID = cpg_positions$Manifest_ID,
    Chromosome = cpg_positions$Chromosome,
    CpG_Position = cpg_center + 1,  # Convert to 1-based for reporting
    Window_Start = pmax(1, cpg_center - window_distance + 1),  # +/- distance, ensure positive, convert to 1-based
    Window_End = cpg_center + window_distance + 1,  # +/- distance, convert to 1-based
    Window_Distance = window_distance,
    Total_Window_Size = 2 * window_distance,
    Genome_Build = cpg_positions$Genome_Build,
    Array_Type = cpg_positions$Array_Type,
    stringsAsFactors = FALSE
  )
  
  # Format for common analysis tools
  windows$UCSC_Format <- paste0(windows$Chromosome, ":", 
                                windows$Window_Start, "-", 
                                windows$Window_End)
  
  cat("Windows created successfully!\n")
  
  return(windows)
}

# Example usage ####
# preload manifest for efficiency
manifest <- get_msa_annotation()

# list cpgs
cpg_list <- c("cg25243766",
              "cg17749946",
              "cg11787522",
              "cg00151250",
              "cg13586425",
              "cg21486834")

# get genomic positions of these CpGs
positions <- lookup_multiple_cpgs(cpg_list, manifest)

# get windows of 3Kbp around genomic positions
windows <- create_mqtl_windows(positions, 3000)


